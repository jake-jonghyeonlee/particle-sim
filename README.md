# particle-sim

층별 전극 마스크와 적층 두께로 **이물 발생 공정(source step)·이물 크기별** 2.5D 기하 접촉 위험
hotspot과 Critical Area(CA)를 계산하는 시뮬레이터다. 명목 조건, open 연결성, 공정 산포,
Fab 가중과 외부 기존 sim 픽셀 비교를 제공한다.

> 출력값은 **기하 접촉 위험(geometric contact risk)**이며 실제 불량확률이 아니다.
> 확률값은 입력한 산포 모델에 조건부이며, 산포가 실제 공정값이 아니면 결과도 실제 Fab
> 불량확률이 아니다.

- 전역 규칙: [SPEC.md](SPEC.md)
- 모듈 구조와 데이터 흐름: [ARCHITECTURE.md](ARCHITECTURE.md)
- 모듈별 계약: [docs/arch/](docs/arch/)
- 2nd 복합 형상 시험 계획과 YAML 입출력 가이드:
  [docs/second-test-and-io-guide.md](docs/second-test-and-io-guide.md)
- 코드 보완 실행 순서와 완료 기준:
  [docs/CODE-IMPLEMENTATION-PLAN.md](docs/CODE-IMPLEMENTATION-PLAN.md)

## 설치

```bash
uv sync                      # Python 3.11+, numpy/scipy/pydantic/zarr/pyarrow/matplotlib 등
uv run pytest                # 단위·통합 테스트
uv run ruff check src tests
```

## 빠른 실행 (합성 데이터)

```bash
uv run particle-sim make-synthetic --case pilot --output work/case
uv run particle-sim validate-input --product-config work/case/manifest.yaml
uv run particle-sim simulate-nominal --product-config work/case/manifest.yaml --output work/run
uv run particle-sim make-report --run work/run --product-config work/case/manifest.yaml

# 공정 산포를 반영한 확률 hotspot (Sobol 2^m 조건)
uv run particle-sim simulate-probability --product-config work/case/manifest.yaml \
    --scenario-config configs/scenarios/nominal.yaml --sobol-m 7 --seed 42 \
    --additional-seed 43 --additional-seed 44 --output work/prob
uv run particle-sim make-report --run work/prob --product-config work/case/manifest.yaml

# 합성 Fab 표로 이전 호환 가중 경로 확인(make-synthetic이 표를 함께 만든다)
uv run particle-sim apply-fab --ca work/prob/ca.parquet \
    --fab-data work/case/fab_particles.csv --product-config work/case/manifest.yaml \
    --bin-method piecewise_constant --mode actual --output work/fab_risk.parquet

# 외부 기존 sim binary hotspot과 픽셀/CA 비교
uv run particle-sim compare-sim \
    --reference-config work/legacy-reference.yaml --run work/run \
    --output work/legacy-comparison
```

합성 case는 `parallel_lines`, `vertical_aligned`, `diagonal_offset`, `triple_overlap`, `pilot`이다.
물리 시뮬레이션 결과가 아니라 코드 검증용 도형이다.

## 입력

제품 manifest YAML 하나가 모든 입력을 연결한다(`docs/arch/schema-config.md`).

- `coordinate_system`: 픽셀 실제 크기(µm/px)와 ROI 크기
- `stack`: 공정 순서, 재질, 두께, 평탄화 계수 — z는 누적으로 계산한다
- `layer_masks`: 레이어 원본 패턴 마스크(전극 분리 전)
- `electrodes`: 레이어별로 분리된 전극 마스크와 net ID. 각 전극은 layer mask의 부분집합이며
  합집합 밖 TR·비평가 영역은 허용하고 면적을 기록한다
- `particle_source_steps`: 이물이 어느 공정 표면에 떨어지는지와 표면 모델(flat/conformal/partial)
- `electrode_pairs` / `excluded_electrode_pairs`: 계산할 조합과 사유가 있는 제외 조합. schema 2.0은
  다른 net의 미분류 조합을 거부한다
- `open_paths`: open 판정 대상 전극, 두 단자 영역 mask와 적용할 이물 유입 공정
- `embedding_depth_um`, `contact_margin_um`, `post_step_policy`: 미확정 물리값. 기본값 없이 명시해야 한다

산포 시나리오는 별도 YAML이다(`configs/scenarios/*.yaml`). `evidence: measured|assumed`가
필수이며 `assumed`면 결과 경고에 그대로 남는다. 모든 산포가 0이면 확률 결과는 nominal과 일치한다.
실측 입력 형식은 `configs/scenarios/fab-overlay-cd.example.yaml`, layer별 평균 이물 표 형식은
`configs/fab-mean-particles.example.csv`를 참고한다. 두 파일의 값은 형식 설명용이다.

마스크는 0/1 또는 0/255 이진 PNG·TIFF만 허용한다. resize·anti-alias된 표시용 이미지는 거부한다.

## 출력

`simulate-nominal` 결과 디렉터리:

| 파일 | 내용 |
|---|---|
| `hotspots.zarr` | 불량 단위별 legacy/nominal hotspot·margin과 mode별 union·동시 단위 수 |
| `ca.parquet` | 전극 조합/open별 raw·`1/k` 배분·고유 CA와 실제 전극/layer ID |
| `ca_union.parquet` | 불량 종류별·전체(`any`) 명목 합집합 CA와 단위별 단순합 |
| `legacy_vs_edt.csv` | 기존 Erode/Dilate와 EDT 기준 hotspot의 pixel 차이와 원인 분류 |
| `run_manifest.json` | config hash, 입력 파일 hash, TODO_PHYSICS 경고, 단위, 버전 |
| `report/` | 그림과 한국어 보고서(`make-report`) |

`simulate-probability` 결과 디렉터리:

| 파일 | 내용 |
|---|---|
| `probability.zarr` | `probability`, `contact_margin_mean_um` `[step, size, unit, y, x]`와 `union_probability` |
| `ca.parquet` | 단위별 nominal CA와 확률 CA 평균·p05·p50·p95, scenario_count, convergence_delta |
| `ca_union.parquet` | 명목 합집합과 산포 평균·p05·p50·p95를 분리한 CA |
| `convergence.csv` | seed별 64/128/256 조건 CA와 직전 대비 상대 변화 |
| `run_manifest.json` | 모든 seed, Sobol m, 시나리오 수, 산포 근거, artifact hash |

`export-simcompare`는 [particle-sim-compare](https://github.com/jake-jonghyeonlee/particle-sim-compare)의
adapter 입력 형식(`bundles.csv`, `layouts.csv`, `hotspots.csv`, `ca_records.csv`,
`layer_stack.csv`, `geometry.json`, `provenance.json`, `arrays/*.npy`)으로 내보낸다.
Fab 자료(`fab.json` / `fab_measurements.csv`)는 생성하지 않으며 `export_warnings.json`에 남는다.

## 핵심 규칙

- 이물 발생 공정과 이물 크기 축을 끝까지 유지한다. 크기를 먼저 합친 hotspot은 원본으로 저장하지 않는다.
- 전체 위험은 pair별 CA 합이 아니라 같은 조건에서 hotspot을 OR한 합집합으로 계산한다.
- 이물보다 나중에 형성되는 layer는 `post_step_policy`로 처리한다.
  `exclude`(엄격) / `protrusion`(돌출·피복 모델) / `wrap_contact`(보수).
- vertical 결과는 보정 전까지 `potential_vertical_short`로만 표기한다.
- 미확정 물리값은 설정 필수값이며 `TODO_PHYSICS` 경고로 결과에 남는다.
- Open은 단순 겹침이 아니라 단자 간 연결성으로 판정한다. 후보 축소는 전수 계산과 동일한 결과를
  내는지 테스트로 확인한다.
- 수렴 허용 오차는 코드에 고정하지 않는다. `convergence.csv`를 보고 사람이 정한다.
- Fab 분포는 geometry 계산이 끝난 뒤에만 적용한다. 크기 bin 정의·단위·조사 면적이 불명확하면
  가중 계산을 중단한다.
