# particle-sim

층별 전극 마스크와 적층 두께로 **이물 발생 공정(source step)·이물 크기별** 2.5D 기하 접촉 위험
hotspot과 Critical Area(CA)를 계산하는 시뮬레이터다. 1차 범위는 Phase 0~3이다.

> 출력값은 **기하 접촉 위험(geometric contact risk)**이며 실제 불량확률이 아니다.
> 산포(Sobol), open 연결성 판정, Fab 분포 가중은 2차 범위다.

- 전역 규칙: [SPEC.md](SPEC.md)
- 모듈 구조와 데이터 흐름: [ARCHITECTURE.md](ARCHITECTURE.md)
- 모듈별 계약: [docs/arch/](docs/arch/)

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
uv run particle-sim export-simcompare --run work/run \
    --product-config work/case/manifest.yaml --output work/simcompare
uv run particle-sim make-report --run work/run --product-config work/case/manifest.yaml
```

합성 case는 `parallel_lines`, `vertical_aligned`, `diagonal_offset`, `triple_overlap`, `pilot`이다.
물리 시뮬레이션 결과가 아니라 코드 검증용 도형이다.

## 입력

제품 manifest YAML 하나가 모든 입력을 연결한다(`docs/arch/schema-config.md`).

- `coordinate_system`: 픽셀 실제 크기(µm/px)와 ROI 크기
- `stack`: 공정 순서, 재질, 두께, 평탄화 계수 — z는 누적으로 계산한다
- `layer_masks`: 레이어 원본 패턴 마스크(전극 분리 전)
- `electrodes`: 레이어별로 분리된 전극 마스크와 net ID. 같은 레이어 전극의 합집합은 레이어 마스크와 일치해야 한다
- `particle_source_steps`: 이물이 어느 공정 표면에 떨어지는지와 표면 모델(flat/conformal/partial)
- `electrode_pairs`: lateral_short / vertical_short 판정 대상
- `embedding_depth_um`, `contact_margin_um`, `post_step_policy`: 미확정 물리값. 기본값 없이 명시해야 한다

마스크는 0/1 또는 0/255 이진 PNG·TIFF만 허용한다. resize·anti-alias된 표시용 이미지는 거부한다.

## 출력

`simulate-nominal` 결과 디렉터리:

| 파일 | 내용 |
|---|---|
| `hotspots.zarr` | `legacy_failure`, `nominal_failure`, `pair_margin_um` `[step, size, pair, y, x]`와 `union_failure`, `union_pair_count` `[step, size, union, y, x]` |
| `ca.parquet` | 전극쌍별 CA(µm², mm², ROI 정규화)와 기존 방식 CA |
| `ca_union.parquet` | 불량 종류별·전체(`any`) 합집합 CA와 pair 단순합 |
| `legacy_vs_edt.csv` | 기존 Erode/Dilate와 EDT 기준 hotspot의 pixel 차이와 원인 분류 |
| `run_manifest.json` | config hash, 입력 파일 hash, TODO_PHYSICS 경고, 단위, 버전 |
| `report/` | 그림과 한국어 보고서(`make-report`) |

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
