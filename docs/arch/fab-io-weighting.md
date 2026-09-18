# arch: fab-io-weighting (P2-T2)

Files: `src/particle_sim/io/fab.py`, `src/particle_sim/aggregation/fab_weighting.py`,
tests `tests/unit/test_fab_weighting.py`.
Depends only on `particle_sim.errors.InputContractError`, numpy and pandas.
Do NOT import geometry/, failure/ or edit any other file.

Fab 자료는 geometry 계산이 끝난 뒤에만 적용한다. 크기 bin 정의·단위·조사 면적이 불명확하면
계산을 멈춘다(fail-closed). 결과는 기대 결함 수이지 실제 불량확률이 아니다.

## io/fab.py
```python
FAB_COLUMNS: tuple[str, ...] = (
    "fab_id", "product_family", "source_step", "size_bin_left_um", "size_bin_right_um",
    "size_representative_um", "particle_count", "inspection_area_mm2", "measurement_period",
    "density_per_mm2", "particle_type",
)
BIN_METHODS: tuple[str, ...] = (
    "exact_point", "piecewise_constant", "piecewise_linear", "empirical_samples",
)

@dataclass(frozen=True)
class FabTable:
    frame: pandas.DataFrame   # FAB_COLUMNS order, density_per_mm2 filled in
    bin_method: str
    source: str               # file path as given

def load_fab_table(path: Path, *, bin_method: str) -> FabTable
```
- `.parquet` -> `pandas.read_parquet`, `.csv` -> `pandas.read_csv`; 다른 확장자는
  InputContractError(file=path, field="fab_data").
- 읽기 실패(OSError/ValueError)는 InputContractError로 변환한다.
- `bin_method`가 `BIN_METHODS`에 없거나 `empirical_samples`이면 로드 시점에
  InputContractError(field="bin_method").
- 수치 열에 비어 있거나 유한하지 않은 값이 있으면 거부한다(NaN 전파 금지).
  기본값을 임의로 고르지 않는다(호출자가 반드시 지정).
- 필수 열 누락 -> InputContractError(field=첫 누락 열). 빈 표 -> InputContractError.
- 검증(모두 InputContractError, field에 열 이름):
  1. `size_bin_left_um < size_bin_right_um`, 둘 다 > 0
  2. `size_bin_left_um <= size_representative_um <= size_bin_right_um`
  3. 같은 `(fab_id, product_family, source_step)` 안에서 bin 구간이 겹치지 않음
  4. `particle_count >= 0`, `inspection_area_mm2 > 0`
  5. `measurement_period`가 비어 있지 않음
  6. `density_per_mm2`가 비어 있으면 `particle_count / inspection_area_mm2`로 채우고,
     값이 있으면 그 계산값과 상대오차 1e-6 안에서 일치해야 한다(불일치 -> 오류)
  7. 같은 `(fab_id, product_family, source_step, size_bin_left_um)` 행 중복 -> 오류
- `product_family`가 비어 있으면 `"*"`(모든 제품군 공통)으로 정규화한다.

```python
def fab_density_for(table: FabTable, *, fab_id: str, product_family: str,
                    source_step: str) -> pandas.DataFrame
```
`fab_id`와 `source_step`이 일치하고 `product_family`가 같거나 `"*"`인 행만 돌려준다
(정확히 일치하는 제품군이 있으면 그 행만 쓰고 `"*"`는 무시). 결과가 비면
InputContractError(field="source_step").

## aggregation/fab_weighting.py
```python
WEIGHTED_COLUMNS: tuple[str, ...] = (
    "run_id", "product_id", "fab_id", "product_family", "weighting_mode", "source_step",
    "failure_mode", "pair_id", "size_bin_left_um", "size_bin_right_um",
    "size_representative_um", "particle_size_um", "ca_um2", "ca_mm2", "density_per_mm2",
    "expected_defect_count", "bin_method",
)

def match_sizes(ca: pandas.DataFrame, density: pandas.DataFrame, *, bin_method: str
                ) -> pandas.DataFrame
```
`ca`는 `particle_size_um`과 `ca_*_um2` 열을 가진다. bin_method별 연결 규칙:
- `exact_point`: `particle_size_um == size_representative_um`인 행만 곱한다(정확 비교, 허용 오차 없음).
- `piecewise_constant`: 반개구간 `left <= size < right`로 판정한다(마지막 bin만 오른쪽 경계 포함).
  경계가 맞닿은 bin에서도 입력 행 순서와 무관하게 같은 결과가 나와야 한다.
- `piecewise_linear`: 대표점 사이 선형 보간으로 각 `particle_size_um`의 density를 구한다.
  관측 범위 밖은 외삽하지 않고 제외한다. 보간에 사용한 양쪽 대표 bin을 모두 연결된 것으로 본다.
- `empirical_samples`: 이 단계에서 지원하지 않는다 -> InputContractError(field="bin_method").
연결되지 않은 (size, bin)은 버리지 말고 반환 프레임의 `unmatched` 속성 대신 별도 반환값으로 준다:
`match_sizes`는 `(matched, unmatched_sizes: list[float], unmatched_bins: list[float])`
튜플을 돌려준다.

```python
def expected_defect_count(ca_mm2: np.ndarray, density_per_mm2: np.ndarray) -> float
# float(np.sum(ca_mm2 * density_per_mm2)); 길이 불일치 -> ValueError
def apply_fab_weighting(ca: pandas.DataFrame, table: FabTable, *, weighting_mode: str,
                        fab_id: str, product_family: str, run_id: str, product_id: str,
                        require_full_coverage: bool = True) -> pandas.DataFrame
```
- `weighting_mode`는 `"common"`(모든 제품에 같은 기준 분포) 또는 `"actual"`(제품에 연결된 Fab).
  다른 값 -> InputContractError(field="weighting_mode").
- `ca`의 각 `source_step`에 대해 `fab_density_for`로 분포를 찾고 `match_sizes`로 연결한 뒤
  `ca_mm2 * density_per_mm2`를 `expected_defect_count` 열에 넣는다.
  `ca_mm2`는 `ca_um2 * 1e-6`로 계산한다.
- `require_full_coverage=True`인데 연결되지 않은 size나 bin이 있으면
  InputContractError(field="particle_size_um")로 중단한다(누락 목록을 메시지에 포함).
- 반환 프레임의 열은 정확히 `WEIGHTED_COLUMNS`, 행 순서는 입력 `ca` 순서를 유지한다.
- geometry 결과(`ca`)는 변경하지 않는다.

Tests (tests/unit/test_fab_weighting.py) — 모두 손계산 가능한 작은 표를 쓴다:
- density 자동 계산(개수/면적)과 기존 값 불일치 거부
- bin 겹침, 대표값 범위 이탈, 음수 개수, 면적 0, 기간 누락, 중복 행, 열 누락 거부
- `bin_method`를 지정하지 않거나 잘못 지정하면 오류, `empirical_samples`는 미지원 오류
- `exact_point`/`piecewise_constant`/`piecewise_linear` 각각의 연결 결과와
  `expected_defect_count` 값이 수작업 계산과 정확히 일치(예: CA 2000 µm² = 0.002 mm²,
  density 5/mm² -> 0.01)
- `require_full_coverage=True`에서 누락 bin이 있으면 중단, `False`면 연결된 것만 계산
- 같은 CA에 common/actual 두 분포를 적용하면 결과가 달라지고 입력 CA 프레임은 그대로다
