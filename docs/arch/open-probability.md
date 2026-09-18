# arch: open, probability and fab integration (phase 2, implemented directly)

이 문서는 Phase 2에서 직접 구현한 모듈의 계약이다. Qwen이 만든 두 leaf 모듈 계약은
`uncertainty-sampler.md`, `fab-io-weighting.md`에 있다.

## failure/open_failure.py

```python
CONFIRMED = "open_confirmed"; PROXY = "open_proxy"

@dataclass(frozen=True)
class OpenResult:
    path_id, electrode, failure (bool y,x), classification, terminal_a_yx, terminal_b_yx, candidates

def blocked_mask(distance_to_center_um, radius_um, contact_margin_um) -> np.ndarray
def causes_open(electrode_mask, blocked, terminal_a_yx, terminal_b_yx) -> bool
def resolve_terminals(open_path, electrode_mask, product_id) -> (terminal_a, terminal_b, classification)
def open_candidates(electrode_mask, grid, radius_um, contact_margin_um, terminals=()) -> np.ndarray
def simulate_open(manifest, masks, open_path, diameter_um, *, contact_margin_um=None) -> OpenResult
```

규칙과 근거:
- open은 겹침이 아니라 **단자 간 연결성**으로 판정한다. 단자가 막힌 경우(label 0)도 open이다.
- 후보는 **전극(또는 단자)에서 도달 거리 안에 있는 중심 전체**다. 이물이 전극을 전혀 덮지 않으면
  연결성을 바꿀 수 없으므로 이 집합은 증명 가능한 superset이다. skeleton과 국부 반폭으로
  줄이는 방식을 먼저 시도했으나 모서리·분기에서 `skeletonize`가 연결을 유지하는 픽셀을 지워
  전수 계산 결과를 놓쳤다(Codex 검토 지적). 등방·비등방 격자, 직선과 ㄱ자 배선, 0이 아닌
  contact margin에서 누락이 없음을 테스트가 확인한다.
- 연결성은 `scipy.ndimage.label` 기본 4-연결을 쓴다(작업지시서 13.1). 대각선 연결을 끊긴 것으로
  보므로 open을 과대평가할 수 있으며, 이는 문서화된 가정이다.
- 단자가 없으면 `fallback`에 따라 `roi_boundary` proxy를 쓰고 결과를 `open_proxy`로 표시한다.
  `skeleton_cut`은 아직 지원하지 않는다(오류).

## units.py
Open 결과를 `PairResult`로 감싸 pair와 동일한 축으로 다룬다. `failure_mode="lateral_open"`,
`electrode_b="N/A"`, margin은 NaN이다. CA 표·union·저장소·simcompare export가 모두 같은 경로를 쓴다.

## geometry/transforms.py
```python
def signed_distance_um(mask, grid)   # 픽셀 중심에서 mask 경계까지 거리(내부 음수). <=0 이 원 mask
def transform_mask(mask, grid, *, overlay_xy_um=(0,0), cd_bias_um=0.0)
def apply_scenario_to_masks(manifest, masks, scenario) -> MaskSet
def apply_scenario_to_stack(stack_z, scenario) -> dict[str, LayerZ]  # 두께 delta와 평탄화 계수
def scenario_manifest(manifest, scenario) -> ProductManifest
```
- overlay는 nearest-neighbour 반복 이동 대신 거리장을 보간해 재임계화한다(작업지시서 7.3).
- CD bias 양수는 확장, 음수는 축소다. 경계는 픽셀 중심에서 반 픽셀 떨어져 있는데, 비등방
  격자에서는 축마다 그 거리가 다르다. 그래서 가장 가까운 반대 위상 픽셀의 방향을 구해
  그 방향의 반 픽셀(`0.5*(|uy|*py + |ux|*px)`)을 뺀다. 임계값 0은 항상 원 mask와 같다.
- 이물 중심 offset은 모든 mask를 반대 방향으로 옮기는 것과 같게 처리한다.

## aggregation/probability.py
```python
class ProbabilityAccumulator: add(failure, margin_um=None); probability(); mean_margin_um()
@dataclass(frozen=True) class ProbabilityCA: ca_probability_mean_um2, ca_scenario_mean_um2,
    ca_p05_um2, ca_p50_um2, ca_p95_um2, scenario_count
def probability_ca(accumulator) -> ProbabilityCA
```
확률은 `failure_count / scenario_count`(float32)이고 CA는 확률 map의 공간 적분이다.
백분위는 scenario별 CA 분포에서 계산하며 단일 평균 map에서 임의로 만들지 않는다.

## validation/convergence.py
`convergence_series({scenario_count: ca})`가 직전 대비 상대 변화를 만든다. 허용 오차는
`is_converged(steps, tolerance)`의 인자이며 코드에 기본값을 두지 않는다.

## pipeline.run_probability
1. manifest + scenario 설정을 읽고 출력 디렉터리가 비어 있는지 확인한다.
2. `build_scenarios`로 Sobol 조건을 만든다(산포 0이면 1개).
3. 조건마다 mask/stack/스칼라를 변환해 `evaluate_units`를 돌리고, 단위별·합집합별로 누적한다.
4. 같은 실행에서 변환 없는 nominal 결과도 계산해 비교용으로 CA 표에 남긴다.
5. Sobol 접두 성질을 이용해 64/128/256 시점의 CA를 `convergence.csv`에 기록한다. 시나리오 수가
   접두에 없으면 전체 수도 함께 기록한다. 단위별과 합집합(`union:<mode>`) 계열을 모두 남긴다.
   `--convergence-tolerance`를 주면 만족하지 못한 계열과 판정 불가 계열을 run manifest에 남긴다.
6. run manifest에 seed, Sobol m, 시나리오 수, 산포 근거 경고를 남긴다.
