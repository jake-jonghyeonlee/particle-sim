# ARCHITECTURE: particle-sim

## 주요 모듈

| 영역 | 모듈 | 역할 |
|---|---|---|
| 입력 | `schema`, `config`, `io/masks`, `io/stack` | schema 1/2 검증, mask와 적층 로드, 잔여 layer 영역 기록 |
| 조합 감사 | `validation/electrode_pairs` | 다른 net 조합을 declared/excluded/unresolved로 분류 |
| 기하 | `geometry/surfaces`, `particles`, `contact`, `transforms` | 표면·구형 이물·3D 거리·overlay/CD 변환 |
| 판정 | `failure/lateral_short`, `vertical_short`, `open_failure`, `units` | short와 단자 영역 open을 공통 불량 단위로 평가 |
| 집계 | `failure/union`, `aggregation/critical_area`, `probability`, `risk_allocation` | OR union, CA, 산포 통계, 1/k 배분 |
| 산포 | `uncertainty/sampler`, `validation/convergence` | 축별 독립 Sobol 조건, 여러 seed, 수렴 기록 |
| Fab | `io/fab`, `aggregation/fab_weighting`, `risk_tables` | layer별 평균 이물 밀도 연결과 네 가지 위험도 관점 |
| 기존 sim | `comparison/legacy_sim` | 외부 binary hotspot 좌표 계약, XOR/CA 비교 |
| 출력 | `io/outputs`, `visualization` | Zarr, Parquet/CSV, manifest, 그림과 실행 보고서 |
| 호환 | `export/simcompare` | 후속 particle-sim-compare용 기존 bundle exporter |

## 명목 실행 흐름

```text
manifest.yaml
  -> schema 검증 + 전극 조합 완전성 감사
  -> layer/electrode/terminal masks + stack z
  -> source step × particle size
       -> 전극별 XY 거리와 Z 거리로 접촉 지도
       -> lateral_short / vertical_short / open 불량 단위 지도
       -> mode별 OR union + 동시 불량 단위 수
       -> raw / 1-k allocated / unique CA
  -> hotspots.zarr + ca.parquet + ca_union.parquet
  -> run_manifest.json (입력·출력 hash, 단위, 가정, 코드 상태)
```

외부 기존 sim 비교는 명목 실행 뒤 독립 단계다. reference 설정이 좌표축 `yx`, 원점 `upper_left`,
픽셀 크기, 직경 의미와 각 binary map의 조건을 명시한다. 공통 불량 단위만 직접 비교하고 기존에 없던
단위는 기준값을 0으로 만들지 않고 `new_pair`로 남긴다.

## 산포와 Fab 실행 흐름

```text
manifest + scenario.yaml
  -> seed별 Sobol 조건 생성
  -> 조건마다 mask/stack 변환 후 불량 단위 평가
  -> 조건 안에서 OR와 1/k 배분
  -> probability.zarr + CA 평균/분위수 + convergence.csv
  -> Fab layer × size 평균 이물 밀도 연결
  -> layer / pair / open / electrode 위험도 + risk_manifest.json
```

Sobol은 균일 공간의 점을 각 입력 분포의 역누적분포로 변환한다. 따라서 NormalSpec 결과는 유한 표본
오차 안에서 정규분포를 따르며, 현재 변수 간 상관은 넣지 않는다. 산포가 가정값이면 결과도 해당 가정에
조건부인 기하 접촉 위험이다.

## 메모리와 저장

2-D 지도 묶음은 chunked Zarr, 분석표는 Parquet, 사람이 볼 순위는 CSV다. 산포 mask가 작으면 메모리
cache를 사용하고 예상 크기가 256 MB를 넘으면 scenario별 bool mask를 packbits 임시 파일에 저장한다.
이 선택과 예상 eager 메모리는 실행 manifest에 기록된다.

세부 입력·출력 예시는 `docs/second-test-and-io-guide.md`, 향후 물리·기능 항목은
`docs/CODE-CHANGE-BACKLOG.md`를 따른다.
