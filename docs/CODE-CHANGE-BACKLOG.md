# particle-sim 코드 수정 백로그

> 실제 구현 순서와 출력 계약은
> [CODE-IMPLEMENTATION-PLAN.md](CODE-IMPLEMENTATION-PLAN.md)를 기준으로 한다. 이 백로그는 각 항목의
> 발견 근거와 상세 설명을 보존한다. 2026-09-19 재검토에서 코드 보완을 먼저 완료한 뒤 기존 sim을
> 비교하는 순서로 확정했다.

작성일: 2026-09-18  
목적: 합성 검증 이후 발견된 구현 보완 사항의 원인, 수정 이유, 대상 코드와 완료 기준을 보존한다.

**2026-09-19 상태:** P0 구현과 자동 검증은 완료했다. 외부 기존 sim 비교 엔진도 구현했으나 실제
현장 산출물이 없어 동일 입력 대조는 아직 실행하지 못했다. P1 simcompare 계약, open 고속화,
tile+halo와 추가 물리는 대체 판정 뒤 진행한다. 아래의 “현재 동작”은 발견 당시 동작을 설명하는
기록이며 현재 코드는 `CODE-IMPLEMENTATION-PLAN.md`의 구현 상태를 따른다.

## 우선순위와 구현 순서

| 순서 | 우선순위 | 작업 묶음 | 적용 시점 |
|---:|---|---|---|
| 1 | P0 | 입력·출력 schema와 현재 회귀 기준 고정 | 모든 코드 변경 전 |
| 2 | P0 | layer 부분집합, 전극 조합 전수 검사, open 단자/적용 공정 | 실제 제품 mask 입력 전 |
| 3 | P0 | overlay·CD 산포, union CA, 중복 없는 기여도와 메모리 보완 | 실제 제품 산포 계산 전 |
| 4 | P0 | Fab 평균 이물 자료와 layer·전극 위험도 출력 | Fab 연계 전 |
| 5 | P0 | 실행 이력·출력 hash와 schema 2.0 고정 | 결과 비교 전 |
| 6 | P0 | 기존 현장 `sim` 비교기와 실제 동일 입력 비교 | 대체 판정 전 |
| 7 | P1 | 실제 제품 pilot와 simcompare 계약 반영 | 기존 sim 비교 후 |
| 8 | P1 | open 고속화와 대형 ROI 분할 계산 | 정확도 동등성 확인 후 |
| 9 | P2 | 결과 그림과 표시 용어 정리, 보고서 재생성 | 결과 고정 후 |
| 10 | 보류 | TR 전기 기능, 재질별 접촉, 비원형 이물 등 추가 물리 | 측정 근거 확보 후 |
| 11 | P0 · 구현 완료 | open 기준 연결성·단자 소실 결함 수정, contact 기반 자동 단자 open, CNO 추가 (12·13절) | 실제 제품 mask 입력 전, 기존 sim 자료 확보와 병행 |
| 12 | P1 · 구현 완료 | 2026-09-19 검토의 비교 계약·수렴 지표·Fab 경로·임시 폴더 결함 (14절) | 기존 sim 대체 판정 전 |

아래 절 번호는 발견 항목을 참조하기 위한 번호이며 구현 순서를 뜻하지 않는다. P0 코드 보완과 전체
회귀 검증을 먼저 끝내고 출력 schema를 고정한 뒤 기존 `sim`과 동일 입력을 비교한다. open 고속화,
대형 ROI 분할과 추가 물리 모델은 대체 비교 이후 별도 작업으로 유지한다.

## 1. 기존 현장 sim 결과 비교기

### 현재 동작

- 저장소의 `legacy_failure`와 대각거리 적용 결과는 비교할 수 있지만, `legacy_failure`는 xy kernel
  차이 확인용 내부 참고값이며 외부 기존 sim 전체를 재현한 결과로 확인되지 않았다.
- 외부 현장 `sim`의 출력 형식 변환기와 픽셀 비교 명령은 없다.
- 따라서 합성 테스트 통과만으로 현장 `sim`을 대체할 수 있는지 판정할 수 없다.

### 수정 이유

좌표 원점, x/y 순서, 픽셀 크기, 직경/반경 해석, 경계 padding과 kernel 모양 중 하나만 달라도
hotspot과 CA가 달라진다. 외부 결과를 같은 축으로 변환해 픽셀 단위로 비교하지 않으면 차이가
2.5D 개선 때문인지 입력 변환 오류 때문인지 구분할 수 없다.

### 수정 대상

- 신규 모듈: `src/particle_sim/comparison/legacy_sim.py`
- CLI 연결: `src/particle_sim/cli.py`
- 입력 계약: 기존 `sim`의 배열 형상, 축, 단위, 전극 조합 ID와 이물 크기 대응표
- 출력: 전극 조합별 XOR 픽셀, 기존 전용/현재 전용 픽셀, CA 절대·상대 차이, 겹침 그림
- 기존 sim과 현재 구현이 공통으로 계산하는 전극 조합만 CA 증감 비교
- G–A와 같이 기존 계산 대상에 없던 대각 전극 조합은 기존값을 0으로 두지 않고 `N/A`와 `new_pair`로 표시
- 비교한 실행의 `run_manifest.json`에 git commit, working tree 변경 여부, 출력 schema version과
  비교기 version을 기록한다. git 정보가 없는 배포 환경에서는 package/build ID를 기록한다.
- 주요 Zarr·Parquet·CSV 출력의 상대 경로와 hash를 실행 또는 비교 manifest에 기록해 전달 후 파일
  교체·누락을 확인할 수 있게 한다.

### 완료 기준

- 동일한 합성 입력에서는 외부 형식으로 내보냈다가 다시 읽어도 XOR가 0이다.
- 축 반전, 단위 불일치, 누락된 전극 조합과 이물 크기를 조용히 보정하지 않고 입력 오류로 낸다.
- 실제 동일 입력에서 차이를 좌표계, 경계 처리, kernel 양자화, kernel 모양, 2.5D 변화로 분류한다.
- 기존 sim에 없는 대각 전극 조합은 새로 계산된 CA로 분리하고 전체 CA 증가는 중복 제거 결과로 확인한다.
- 미분류 차이가 남으면 대체 판정을 실패로 기록한다.
- git 저장소에서 실행한 결과의 `git_commit`이 `null`로 남지 않고, 수정 파일이 있으면 dirty 상태가
  함께 기록되어 같은 코드 버전을 식별할 수 있다.

## 2. layer mask와 전극 영역의 부분집합 관계 허용

### 현재 동작

`src/particle_sim/io/masks.py`는 한 layer의 모든 전극 mask 합집합이 layer mask와 정확히 같지 않으면
입력을 거부한다. layer mask에 있지만 전극 mask에는 없는 픽셀도 `missing`으로 처리한다.

### 수정 이유

실제 layer mask에는 평가 전극 외에 TR 역할 영역이나 기타 패턴이 포함될 수 있다. 따라서 전극 mask
합집합은 layer mask의 일부일 수 있으며 완전 일치를 강제하면 정상적인 실제 제품 입력을 읽지 못한다.

TR은 전기 신호를 분리하는 역할을 하지만, 현재 이물 해석은 그 전기적 기능을 모델링하지 않는다.
따라서 이번 변경에서는 남는 TR 영역을 허용하고 그 존재와 면적만 기록한다. TR이 이물에 의해 영향을
받았을 때 신호 분리 기능이 어떻게 변하는지, 그리고 이를 short·open 또는 별도 불량 유형으로 판정할지는
물리 기준과 검증 데이터가 정해진 뒤 별도 개발 항목으로 결정한다.

### 수정 대상

- `src/particle_sim/io/masks.py`
  - 모든 전극 mask가 대응하는 layer mask의 부분집합인지 검사
  - 같은 layer의 전극 mask끼리 중첩되지 않는 기존 검사는 유지
  - `layer mask - 전극 mask 합집합` 영역은 허용하고 픽셀 수·비율을 검증 결과에 기록
- `src/particle_sim/schema.py`
  - 남는 영역을 TR/기타 미평가 영역으로 표시할 수 있는 선택 메타데이터 검토
  - 이 메타데이터는 현 단계에서 불량 판정에 사용하지 않고 추후 모델 확장 근거로만 보존
- `validate-input`
  - 전극이 layer 밖으로 나간 픽셀은 오류
  - layer 안의 미할당 영역은 픽셀 수·면적 비율과 함께 정상 정보로 출력

### 완료 기준

- 전극 mask가 layer mask 밖으로 나가면 입력 오류가 발생한다.
- layer mask에 TR 영역이 남아 있어도 입력을 통과하고 남는 픽셀 수·면적 비율을 기록한다.
- 현재 short·open·hotspot 판정은 전극 mask를 기준으로 하며 TR의 전기적 신호 분리 기능은 반영하지 않는다.
- 결과 manifest와 보고서에 이 제한을 기록하여 추후 TR 모델 검토 시 현재 결과와 구분할 수 있다.
- 기존 합성 fixture는 검증 규칙을 바꾼 뒤에도 동일한 hotspot·CA를 유지한다.

### 추후 검토 항목 — TR의 전기적 역할

- TR이 정상일 때 분리하는 신호와 연결 관계를 입력으로 어떻게 표현할지
- 이물이 TR을 덮거나 가로지를 때 발생하는 물리 현상과 불량 판정 기준
- 기존 `lateral short`, `vertical short`, `open` 중 어느 유형으로 표현할지 또는 별도 유형이 필요한지
- 회로/소자 해석 또는 Fab 불량 데이터와 결합해야 하는지
- TR 효과를 포함했을 때 기존 sim과 비교 가능한 기준 사례가 있는지

단계별 제안(`tr_hit` 기하 기록 → S–D short·gate short·채널 open 규칙 → kill ratio 보정)과 참고 문헌은
[`report/REVIEW-PLAN.md`](report/REVIEW-PLAN.md)의 "PDF 제작 시 추가할 향후 개선 항목"에 정리했다.

## 3. 전극 조합 누락 방지

### 현재 동작

`src/particle_sim/schema.py`는 `electrode_pairs`에 선언한 전극 조합의 유효성은 검사하지만 가능한
조합을 자동으로 만들거나 선언 누락을 알려주지 않는다. 한 layer에 여러 전극이 있으면 사용자가
목록에서 빠뜨린 조합은 계산 결과에도 나타나지 않는다.

### 수정 이유

복잡한 실제 형상에서는 육안으로 인접 전극 조합을 고르면 lateral short 또는 vertical short 후보를
누락할 수 있다. 누락된 조합은 hotspot이 0인 경우와 출력상 구분되지 않아 구조 위험을 과소평가한다.

### 수정 대상

- 신규 후보 생성/검사 모듈: 예) `src/particle_sim/validation/electrode_pairs.py`
- `validate-input` 출력: 선언 조합, 기하학적 도달 가능 후보, 제외 조합과 제외 사유
- YAML 선택 항목: 수동 승인 목록 또는 후보 생성 규칙과 최대 검토 거리
- 동일 net 제외, 같은 layer는 lateral short, 다른 layer는 vertical short라는 기존 규칙 유지

### 완료 기준

- 여러 전극이 있는 각 layer와 layer 사이에서 후보 조합을 결정적으로 생성한다.
- 선언 목록과 후보 목록의 차이를 오류 또는 명시적 승인 경고로 보여준다.
- 의도적인 먼 거리 음성 대조와 사용자가 승인한 제외 목록을 보존한다.
- 기존 manifest는 명시적인 호환 모드에서 같은 계산 결과를 낸다.

### open 단자 입력 개선

현재 `open_paths`는 연결을 확인할 두 단자를 각각 한 개의 `(y, x)` 픽셀 좌표로 입력한다. 복잡한
실제 배선에서 대표 픽셀 하나를 사용자가 직접 찾기 어렵고, mask 해상도나 형상이 바뀌면 좌표가 전극
밖으로 나갈 수 있다. 단자는 open의 의미를 정하는 데 필요하지만 점 좌표일 필요는 없다.

실제 제품 입력은 두 개의 **단자 영역 mask**를 권장 방식으로 추가한다.

- `terminal_a_mask`, `terminal_b_mask`: 배선 mask와 같은 shape의 작은 이진 영역
- 각 단자 영역은 비어 있지 않고 대상 전극의 부분집합이며 서로 겹치지 않아야 한다.
- 이물이 단자 영역 전체를 가리거나 손상 후 두 단자 영역 사이에 연결 경로가 없으면 open이다.
- 이 규칙은 사용자가 직접 지정한 `open_paths` 단자 영역에 적용한다. contact·via에서 자동으로 만든
  단자가 이물에 막히는 경우는 open이 아니라 CNO로 분류한다(13절).
- 기존 단일 픽셀 좌표는 합성 시험과 이전 manifest 호환용으로 유지한다.
- `roi_boundary` 자동 선택은 실제 제품 결과가 아닌 `open proxy`로만 허용한다.
- 입력 검증 그림에 배선과 두 단자 영역을 함께 표시해 좌표 숫자를 보지 않고도 확인할 수 있게 한다.

완료 기준은 직선·우회·분기·ㄱ자 배선에서 단일 픽셀 단자와 그 픽셀을 포함한 단자 영역 mask가 같은
결과를 내고, 영역의 일부가 가려진 경우와 전체가 가려진 경우의 규칙이 전수 계산과 일치하는 것이다.

## 4. mask layer별 Fab 실측 overlay·CD 산포 입력

### 확인된 입력 자료

mask layer별 overlay와 CD 산포는 Fab 실측값이 있다. 실제 제품 산포 계산에서는 합성 예제의
가정값을 재사용하지 않고 이 측정값을 layer별로 연결한다. 두께, 이물 직경, 이물 중심 offset,
매립 깊이와 접촉 여유는 별도 근거가 필요하므로 overlay·CD와 같은 수준의 실측값으로 취급하지 않는다.
현재 확인 기준으로 mask별 CD와 overlay는 서로 독립이고, overlay x와 y도 독립이다.

### 현재 동작

- `ScenarioConfig.overlay_um`은 layer마다 `NormalSpec` 하나를 받으며, 같은 평균·표준편차를 x와 y에
  각각 독립적으로 적용한다.
- `ScenarioConfig.cd_bias_um`은 layer마다 정규분포 하나를 받는다.
- `evidence`는 시나리오 전체에 `measured` 또는 `assumed` 하나만 기록한다.
- 따라서 overlay·CD는 실측이고 다른 항목은 가정인 혼합 입력을 정확히 표현할 수 없다.
- 실측 자료의 x/y별 산포, 독립 가정의 근거, 기준 layer, lot·공정·기간과 원자료 식별자를 저장할 필드가 없다.

### 수정 이유

overlay는 기준 layer에 대한 상대 이동이며 x와 y의 분포가 다를 수 있다. 독립 관계는 확인됐지만 현재처럼
같은 정규분포를 두 축에 적용하면 Fab 실측 분포와 다른 mask 이동을 만들 수 있다. CD도 실측값의
평균 bias와 산포, 측정 방향과 집계 범위를 잃으면 다른 layer 데이터가 잘못 연결될 수 있다.
시나리오 전체를 `measured`로 표시하면 아직 가정인 두께·이물 항목까지 실측처럼 보이는 문제도 있다.

### 수정 대상

- `src/particle_sim/schema.py`
  - overlay를 layer별 x/y 분포로 각각 표현하고 현재 독립 관계를 메타데이터에 기록
  - 각 산포 항목에 `evidence`, 원자료 ID, Fab 공정, 제품군, lot/기간, 단위와 버전을 기록
  - 기준 layer와 상대 overlay의 부호 규칙을 명시
- `src/particle_sim/uncertainty/sampler.py`
  - overlay x/y와 mask별 CD·overlay를 각각 독립된 Sobol 차원으로 변환
  - Sobol의 0~1 저불일치 표본을 정규분포 또는 실측 경험분포의 역누적분포로 변환
  - 실측 요약 통계와 원본 표본 입력을 구분
- `src/particle_sim/geometry/transforms.py`
  - 같은 mask layer의 layer mask와 모든 전극 mask에 동일한 overlay·CD 변환 적용 유지
  - 기준 layer에는 상대 overlay 0이 적용되는지 검사
- 신규 Fab 산포 입력 변환기
  - 사내 실측 표를 scenario YAML 또는 고정된 중간 형식으로 변환
  - layer 이름 매핑, 단위 변환과 누락값을 검증하고 변환 이력을 보존

### 먼저 확정할 데이터 계약

1. overlay가 절대값인지 기준 layer 대비 상대값인지와 기준 layer
2. x/y 각각의 평균·표준편차; x-y와 mask별 CD·overlay는 독립으로 확정
3. CD가 선폭 전체 변화인지 경계 한쪽 이동량인지와 양수 부호의 의미
4. 실측값의 단위, 제품군, 공정, 장비, lot와 집계 기간
5. 정규분포 적합값을 사용할지 원본 표본을 직접 표본화할지

### 완료 기준

- mask layer별 실측 overlay·CD가 지정한 layer에만 적용된다.
- 같은 layer의 layer mask와 전극별 mask가 같은 변환을 받아 합집합 관계를 유지한다.
- x/y 산포가 다른 시험에서 각 축의 표본 통계가 입력값과 허용 오차 안에서 일치하고 독립 차원으로 처리된다.
- overlay·CD는 `measured`, 아직 근거가 없는 항목은 `assumed`로 각각 결과 manifest에 기록된다.
- 실측 입력 파일 ID·hash·단위·매핑 버전으로 같은 산포 조건을 재현할 수 있다.

## 5. 산포 결과의 중복 제거 전체 CA 열 분리

### 현재 동작

`src/particle_sim/pipeline.py`의 확률 계산은 산포가 반영된 중복 제거 평균 CA를
`ca_union_nominal_um2`에 기록한다. `src/particle_sim/io/outputs.py`의 열 이름은 기준 조건 결과처럼
보이므로 실행 종류를 함께 보지 않으면 의미를 잘못 해석할 수 있다.

### 수정 이유

기준 조건 CA와 산포 평균 CA는 계산 정의가 다르다. 같은 열 이름을 재사용하면 Fab 가중, 보고서,
simcompare에서 기준 조건 값과 산포 평균값을 섞을 수 있다.

### 수정 대상

- `src/particle_sim/io/outputs.py`: 확률용 열 추가
  - `ca_union_mean_um2`
  - `ca_union_p05_um2`, `ca_union_p50_um2`, `ca_union_p95_um2`
  - `scenario_count`, `convergence_delta`
- `src/particle_sim/pipeline.py`: 기준 조건과 산포 결과를 각각 맞는 열에 기록
- `src/particle_sim/visualization/reports.py`: 실행 종류에 따라 올바른 열 표시
- schema version을 올리고 이전 결과를 읽을 때의 호환 규칙 명시

### 완료 기준

- 산포 0에서 `ca_union_mean_um2 == ca_union_nominal_um2`가 성립한다.
- 산포가 있는 실행에서 기준 조건 CA와 산포 평균 CA가 서로 다른 열에 함께 보존된다.
- 중복 제거는 각 산포 조건에서 전극 조합별 hotspot을 먼저 OR한 뒤 평균한다.
- 이전 schema 결과는 명시적인 변환 또는 이해하기 쉬운 오류로 처리한다.

## 6. Fab 공정 연결과 중복 제거 전체 CA 가중

### 현재 동작

- `src/particle_sim/aggregation/fab_weighting.py`는 `ca.parquet`의 `pair_id`를 필수로 요구하고
  전극 조합별 CA 행을 가중한다.
- `source_step` 값으로 Fab 표를 직접 조회하며 manifest의 `fab_distribution_layer`를 사용하지 않는다.
- 전극 조합별 기대값을 합하면 같은 이물 중심 위치가 여러 전극 조합을 만들 때 중복 집계된다.
- 한 Fab 크기 구간에 여러 시뮬레이션 직경을 연결하면 구간 전체 밀도가 반복 적용될 가능성이 있다.
- 현재 `ca_union.parquet`의 `failure_mode=any`는 특정 전극 layer의 전극 형상 합집합이 아니다. 같은
  이물 유입 공정·이물 크기에서 선언된 모든 lateral short·vertical short·open hotspot의 OR다.
- 현재 Fab 입력은 합계 `particle_count / inspection_area_mm2` 형식만 직접 표현한다. 실제 보유 자료는
  layer별·이물 크기 구간별 **평균 이물 개수**이므로 입력 의미가 맞지 않는다.
- 전극 조합별 원본 가중값은 어느 조합이 민감한지는 보여주지만, 같은 위치에서 여러 조합이 동시에
  hotspot이면 합계가 중복된다. 반대로 union 가중값만 남기면 어느 전극 조합이 주요한지 알 수 없다.

### 수정 이유

구조 전체 비교의 대표값은 같은 위치를 한 번만 세는 중복 제거 전체 CA여야 한다. 또한 시뮬레이션의
이물 유입 공정 ID와 Fab 측정 공정 이름은 같다는 보장이 없다. 잘못 연결하면 다른 공정의 이물 밀도를
적용하고도 정상 결과처럼 보일 수 있다.

Fab 자료가 layer별 평균 이물 개수라면 각 layer를 해당 이물 유입 공정에 명시적으로 연결하고, 그
공정의 `failure_mode=any` union CA에 크기별 평균 밀도를 적용한다. 이 union에는 그 공정에서 유입된
이물이 만들 수 있는 이후 layer와의 vertical short도 포함돼야 하므로 “해당 layer에 그려진 전극만의
union”으로 제한하면 안 된다.

최종 목적에는 layer별 전체 위험, 전극 조합별 위험, layer와 전극 조합을 함께 본 위험이 모두 필요하다.
따라서 전체값은 union으로 중복을 제거하고, 같은 union을 전극 조합에 배분한 기여도도 함께 계산해야
한다. 원본 전극 조합 CA는 민감도 순위를 위해 별도로 보존한다.

### 수정 대상

- `src/particle_sim/aggregation/fab_weighting.py`
  - 전극 조합별 `ca.parquet`와 전체 `ca_union.parquet` 입력을 명시적으로 구분
  - 구조 비교는 `failure_mode=any`의 `ca_union_mean_um2`를 기본 입력으로 사용
  - `fab_distribution_layer`로 Fab 표 조회
  - 한 크기 구간의 밀도를 둘 이상의 대표 직경에 중복 적용하지 못하게 검사
- `src/particle_sim/io/fab.py`
  - 조회 인자를 시뮬레이션 ID가 아닌 Fab 측정 layer ID로 명확히 변경
  - layer·크기 구간별 `mean_particle_count`와 `inspection_area_per_sample_mm2`를 입력받아
    `mean_density_per_mm2 = mean_particle_count / inspection_area_per_sample_mm2`로 변환
  - 평균이 이미 면적 정규화된 값이면 별도의 `mean_density_per_mm2` 단위로 입력하고 개수와 혼용 금지
  - 표본 수, 측정 기간과 평균 산출 범위를 추적 정보로 보존
- `src/particle_sim/cli.py`: `weighted_CA_index`와 `expected_defect_count` 모드를 분리
- 결과 계약과 단위 갱신: 상대 비교 지수와 실제 밀도 적용 기대 건수를 다른 열로 저장
- `ca_union.parquet`에 union 범위를 명시한다.
  - 주값: 이물 유입 공정별 모든 선언 불량의 `any` union
  - 원인 분석: lateral short의 전극 layer별 union, vertical short의 layer 조합별 union을 선택 출력
- `risk_by_layer.parquet` 추가
  - `fab_distribution_layer × particle_size`별 `any` union CA, 평균 이물 밀도와 가중 위험
  - 크기 합계와 layer별 순위, 전체 대비 기여율
- `risk_by_pair.parquet` 추가
  - `fab_distribution_layer`, `source_step`, `particle_size`, `pair_id`, 실제 두 전극 ID·layer,
    `failure_mode`를 구조화된 열로 저장
  - `raw_pair_risk`: 전극 조합 hotspot 전체를 인정한 비가산 민감도
  - `allocated_pair_contribution`: 같은 픽셀에서 활성 조합이 `k`개면 각 조합에 `1/k` 배분
  - `unique_pair_contribution`: 해당 조합을 제외했을 때 union이 감소하는 양
  - 산포 계산은 각 scenario에서 배분한 뒤 평균·분위수를 계산
- `risk_by_electrode.parquet` 추가
  - `involvement_risk`: 해당 전극이 포함된 모든 전극 조합의 배분 기여를 모은 비가산 관여 점수
  - `allocated_electrode_contribution`: short는 조합 기여의 절반씩, open은 해당 전극에 전부 배분해
    전극 합계가 union과 일치하는 보고용 기여도
- 사람이 바로 열 수 있는 `risk_summary.csv`에 layer·전극 조합·개별 전극 상위 순위와 단위 출력
- `1/k` 및 전극 절반 배분은 물리적 인과가 아니라 중복 없는 보고 규칙임을 schema metadata에 기록

### 완료 기준

- 겹치는 전극 조합 합성 도형에서 전극 조합별 단순합과 중복 제거 전체 가중값이 다르게 나온다.
- `source_step != fab_distribution_layer`인 시험에서 지정한 Fab 공정의 분포만 사용한다.
- 한 크기 구간을 여러 직경에 연결하면 명시적 배분 규칙이 없는 한 입력 오류가 발생한다.
- 공통 이물 분포를 적용한 구조 순위와 실제 Fab 분포를 적용한 기대값을 별도 결과로 저장한다.
- Fab layer별 평균 이물 개수와 검사 면적을 입력한 수기 계산이 프로그램의 평균 밀도와 일치한다.
- `source_step → fab layer` 연결 후 `union:any × mean_density_per_mm2`의 layer·크기별 기여와 전체 합이
  수기 계산과 일치한다.
- 겹침이 없는 두 전극 조합은 `allocated_pair_contribution == raw_pair_risk`다.
- 완전히 같은 hotspot을 가진 두 전극 조합은 각 배분 기여가 union의 절반이고 두 값의 합은 union이다.
- 부분 중첩 fixture에서 픽셀별 `1/k` 수기 계산과 전극 조합별 배분 기여가 일치한다.
- 산포가 있는 fixture에서 모든 전극 조합의 평균 배분 기여 합이 scenario별 OR 후 평균한 union과 같다.
- short만 있는 fixture에서 개별 전극의 `allocated_electrode_contribution` 합은 union과 같고,
  `involvement_risk`는 합산 불가 지표로 표시된다.

## 7. simcompare 공정 의미 계약 반영

### 현재 동작

`src/particle_sim/export/simcompare.py`는 vertical short일 때 `assigned_process_layer`를 비우고 경고를
남긴다. simcompare 쪽은 상부 전극 layer를 기대하고, particle-sim은 이물 유입 공정과
`fab_distribution_layer`를 기준으로 삼기 때문이다.

또한 particle-sim의 `ca.parquet`에는 `pair_id`, `electrode_a`, `electrode_b`가 있지만 현재 exporter는
이를 simcompare의 구조화된 필드로 전달하지 않는다. `sample_id` 문자열에는 조합 ID가 들어가지만,
simcompare의 비교 그룹과 순위는 `failure_type + 전극 역할 집합 + 전극 layer 집합`을 중심으로 한다.
따라서 동일 역할을 가진 실제 전극이 여러 개이면 “어느 실제 전극 조합이 주요한가”를 안정적으로
집계할 수 없다.

simcompare는 `average_count_per_glass`와 검사면적을 입력 계약으로 읽을 수 있으나, 현재 CAA 해석은
demo mode의 `probability`만 지원한다. 실측 평균 이물 개수를 적용한 전극 조합 순위는 아직 계산할 수
없다. exporter도 현재 `hotspots.zarr`의 기준 조건만 읽고 `probability.zarr` 산포 결과는 내보내지 않는다.

### 수정 이유

빈 값은 잘못된 자동 연결을 막지만 연동을 완료한 것은 아니다. 양쪽 도구가 다른 의미의 공정 값을
같은 필드로 사용하면 내보내기 성공 후에도 잘못된 Fab 분포가 적용될 수 있다.

전극 역할과 layer는 서로 다른 제품의 유사 구조를 비교할 때 필요하고, 실제 전극 ID와 조합 ID는 한
제품 안에서 주요 전극을 찾을 때 필요하다. 두 종류의 식별자를 모두 보존해야 layer별·전극별·교차
위험을 같은 데이터에서 계산할 수 있다.

### 수정 대상

- 양쪽 데이터 계약에서 다음 값을 분리한다.
  - 이물 유입 공정
  - Fab 측정 공정
  - 하부/상부 전극 layer
- simcompare `AnalysisBundle`과 CSV 계약에 `pair_id`, `electrode_a_id`, `electrode_b_id`,
  `particle_source_step`, `fab_distribution_layer`를 추가한다. 전극 역할과 layer 필드는 그대로 유지한다.
- 한 제품 안의 전극 중요도 분석은 `fab_distribution_layer × pair_id × failure_type`을 기본 행으로 하고,
  `allocated_pair_contribution`을 기본 순위값으로 사용한다. `raw_pair_risk`와
  `unique_pair_contribution`도 함께 표시한다.
- layer별 순위, 실제 전극 조합별 순위, 개별 전극별 순위와 `layer × pair` 교차표를 각각 출력한다.
- `average_count_per_glass / surveyed_area`를 크기별 면적 밀도로 해석하는 real-mode 규칙을 명시하고,
  단위와 전체 크기 구간 포함 여부를 검증한다.
- 기준 조건 hotspot과 산포 probability map을 모두 내보내되 `value_semantics`를 구분한다.
- `src/particle_sim/export/simcompare.py`와 simcompare importer를 같은 계약 버전으로 갱신
- 계약 버전과 매핑표를 export manifest에 기록

### 완료 기준

- 이물 유입 공정과 상부 전극 layer가 다른 합성 예제를 양쪽 도구에서 같은 Fab 분포에 연결한다.
- 누락 또는 모호한 매핑은 경고 후 빈 값으로 진행하지 않고 연동 단계에서 실패한다.
- 내보내기→simcompare 읽기→Fab 가중 결과가 수기 계산과 일치한다.
- 동일 역할의 전극이 여러 개인 fixture에서 실제 전극 ID별 순위가 유지되고 `sample_id` 문자열 파싱에
  의존하지 않는다.
- simcompare의 전극 조합별 배분 기여 합과 particle-sim의 layer별 union 가중 위험이 허용 오차 안에서
  일치한다.
- 실측 `average_count_per_glass` 입력으로 real mode 순위가 생성되고 결과가 실제 불량률이나 수율로
  잘못 표시되지 않는다.

## 8. open 계산 캐시와 성능 개선

### 현재 동작

- `src/particle_sim/units.py`는 이물 유입 공정마다 `simulate_open`을 다시 호출한다.
- `src/particle_sim/failure/open_failure.py`는 팽창 연산으로 후보를 만든 뒤 각 후보마다 전체 배선의
  연결 성분을 다시 계산한다.
- 1000×1000 이상, 큰 이물, 많은 크기와 산포 조건에서 후보 수만큼 전체 이미지 라벨링이 반복된다.

### 수정 이유

open은 short와 계산 특성이 다르며 현재 구조에서는 전체 TAT를 지배할 가능성이 높다. 다만 창 단위
계산을 단순 적용하면 창 밖의 우회 경로를 잃어 잘못된 open을 만들 수 있으므로 동등성 증명이 먼저다.

### 수정 대상

- 이물 유입 공정 간 재사용 캐시: 전극 mask, 단자, 이물 크기, 접촉 여유가 같을 때만 재사용
- 후보 생성: 큰 binary dilation 대신 거리변환 기반 후보 검토
- 연결성: 안전한 전역 연결 정보와 국부 변경 영역을 결합하는 방법 검토
- 크기 단조성 이월은 작은 전체 픽셀 기준 계산으로 먼저 증명한 뒤 적용
- benchmark에 open 단독 d=1/10, open 포함 기준 조건, 32크기×64조건 추가

### 완료 기준

- 직선·우회·분기·ㄱ자, 등방·비등방 격자의 모든 이물 중심에서 기존 전수 계산과 픽셀 일치한다.
- 캐시 키가 달라지는 입력에서는 결과를 재사용하지 않는다.
- 같은 장비와 입력에서 수정 전후 TAT와 최대 RSS를 5회 반복 측정해 보존한다.
- 시간 단축과 함께 hotspot, CA, 분류 결과가 완전히 일치한다.

## 9. 대형 ROI 분할 계산

### 현재 동작

1000×1000 전체 이미지는 실행됐지만 더 큰 ROI를 tile+halo로 나누는 계산은 없다. Zarr 출력은
chunk를 사용하더라도 계산 중간 배열은 전체 이미지 크기의 영향을 받는다.

### 수정 이유

실제 입력이 1000×1000보다 크거나 전극 조합·산포 조건이 늘면 메모리가 급격히 증가한다. 단순 타일
분할은 이물 반경과 대각 접촉 거리, open 우회 경로가 타일 경계를 넘을 때 결과를 잘라낼 수 있다.

### 수정 대상과 완료 기준

- short는 최대 이물 반경과 접촉 여유를 포함한 halo 크기를 계산한다.
- 타일 내부 결과를 원래 ROI 좌표로 합치며 경계 중복 픽셀을 한 번만 기록한다.
- open은 전역 연결성이 필요한 별도 알고리즘이 마련되기 전까지 타일 모드에서 거부한다.
- 전체 이미지 계산과 타일 계산의 hotspot·CA가 경계 픽셀까지 일치해야 한다.
- 1000²와 대표 실제 크기에서 최대 RSS 감소를 측정한다.

## 10. 결과 그림과 표의 표시 용어 정리

### 현재 동작과 이유

`src/particle_sim/visualization/hotspots.py`와 `reports.py`가 `pair`, `source step`, `nominal`,
`union`, `legacy`를 그대로 표시한다. 파일·열 이름은 계약상 유지해야 하지만 사람이 읽는 제목과
범례는 전극 조합, 이물 유입 공정, 기준 조건, 중복 제거 전체, 기존 방식으로 표시하는 편이 명확하다.

### 완료 기준

- 화면 표시에서는 `lateral short`, `vertical short`, `open`, `hotspot`을 그대로 사용한다.
- 나머지 내부 용어는 쉬운 한국어로 표시하고, 정확한 파일명과 열 이름은 부록에 함께 적는다.
- 데이터 파일의 schema와 프로그램용 식별자는 이 표시 변경만으로 바꾸지 않는다.

## 11. 요구와 근거 확보 후 결정할 항목

다음은 현재 코드 결함으로 확정하지 않는다. 실제 사용 요구와 물리 근거가 생기면 각각 별도 사양과
검증 기준을 만든 뒤 구현한다.

- `empirical_samples`: 원본 이물 크기 표본을 CA 크기 축에 연결하는 적분 규칙
- `skeleton_cut`: 우회 경로와 분기에서 open 의미를 보존하는 판정 규칙
- vertical short 전기 보정: 전도성, 절연막 균열, 접촉저항과 피복 파괴 근거
- 이물 재질별 전기 접촉 모델: 전도성·비전도성 이물에 따라 short 성립 여부와 유효 접촉거리가
  달라지는지 검토한다. 기하 거리 보정과 전기적 연결 조건을 하나의 `contact_margin_um`에 섞지 않고,
  이물 재질·계면 상태·불량 유형별 판정값으로 분리할 근거를 먼저 확보한다.
- 비원형 이물 형상 모델: 현재 원형 footprint와 구형 높이 모델 외에 타원형, 길쭉한 형상과 불규칙
  형상을 지원할 필요가 있는지 검토한다. 형상별 장·단축, 방향각, 3D 높이와 Fab 형상 분포를 입력으로
  정의하고 lateral short·vertical short·open 결과가 회전 방향에 따라 달라지는지를 검증해야 한다.
- 이물 크기별 Fab 분포: 작은 이물이 많고 큰 이물이 적은 실측 경향을 크기 구간별 개수와 조사 면적으로
  정규화해 CA 계산 후 가중한다. 이 분포를 Sobol의 이물 직경 산포와 섞지 않으며, 공정·위치·기간과의
  상관관계가 확인되기 전에는 공통 분포로 구조 상대 비교에 사용한다.
- 현재 실행의 `contact_margin_um`은 0 µm로 고정한다. 비영값은 실측 또는 물리 근거와 검증 사례가
  마련된 뒤 사용한다.
- 실제 불량확률/수율 변환: 노출 면적·기간, 이물 독립성, Poisson 가정과 제품 대조
- CNO 가림 기준: contact가 **전부** 가려져야 CNO인지, **일부**만 가려져도 저항 증가나 접촉 불량으로
  CNO인지 논의가 필요하다. 1차 구현은 전부 가림을 잠정 기준으로 두고 `TODO_PHYSICS` 경고를 남긴다.
  일부 가림을 채택하면 가림 비율 임계값과 그 근거(접촉저항 측정, FA 사례)를 함께 정한다. 기준을
  재계산 없이 바꿀 수 있도록 contact별 최대 가림 비율을 출력할지도 함께 검토한다.
- CNO 중복 contact: 같은 두 전극을 여러 contact hole이 병렬로 잇는 경우, hole 하나만 막혀도 불량인지
  모두 막혀야 불량인지 정한다. 1차 구현은 hole 하나라도 막히면 CNO로 보는 보수적 기준이다.
- net 단위 open: 같은 net이 다른 layer와 contact를 거쳐 우회 경로를 가지면 전극 하나가 끊겨도 실제
  신호는 유지될 수 있다. 1차 자동 단자 open은 전극 단위로 판정하고, 여러 layer를 묶은 net 연결
  그래프 판정은 후속으로 검토한다.
- 전류 방향과 불량 심각도: contact·via와 적층 순서로는 연결 관계만 알 수 있다. 구동단·부하 정보를
  추가하면 끊긴 쪽이 신호를 잃는지 구분할 수 있으나 open 판정 자체에는 필요하지 않다.

## 12. contact 기반 자동 단자 open

**2026-09-19 구현 완료.** 아래 "현재 동작"은 구현 전 기록이다. 구현 결과와 측정은
`CODE-IMPLEMENTATION-PLAN.md` 단계 5a 상태에 있다. 14절 5번(비교기 `new_unit`·`not_applicable`
분리)도 이 작업에서 함께 반영했다.

### 현재 동작

- open은 `open_paths`에 적은 경로만 계산하며, 단자는 사용자가 영역 mask·점 좌표로 직접 지정한다.
- 입자가 없는 상태에서 두 단자가 이미 연결돼 있는지 확인하지 않는다. 입력 전극이 끊겨 있거나 산포
  변형으로 가는 선이 끊기면 후보 위치 전부가 open으로 계산된다(2026-09-19 검토 1번, 재현: 두 섬 전극
  608/608 후보가 open).
- 산포에서 단자 영역은 overlay만 적용하고 CD는 적용하지 않는다. 단자가 전극 가장자리에 있고 CD가
  줄면 살아남은 단자가 0이 되어 모든 후보가 open이 된다(검토 2번, 재현: 1 px 단자, CD -0.3 µm).

### 수정 이유

실제 layout은 전극 수가 많아 모든 배선의 단자를 사람이 지정하기 어렵다. 단자 없이 연결 성분 수
증가로 판정하는 방식도 검토했지만, 파일럿 G1·D1·D2에서 작은 모서리 조각과 기능 없는 가지 절단을
open으로 세고(직경 10 µm G1 +474곳), 단자가 통째로 덮이는 경우를 놓쳤다(−144곳). 크기 임계값으로는
오검출과 누락이 서로 바뀔 뿐이다. 전극이 기능을 갖는 위치는 다른 layer와 이어지는 contact·via와 ROI
밖으로 이어지는 경계이므로, 이 위치를 자동 단자로 쓰면 입력 부담 없이 open의 의미를 유지할 수 있다.

### 수정 대상

- schema 2.0 `contacts` 입력
  - `contact_id`, `path`(contact·via 원본 mask), `hole_layer`(구멍이 뚫리는 절연층)
  - `lower_layer`·`upper_layer`는 생략하면 stack `process_order`에서 `hole_layer` 바로 아래·위의
    patterned 도체 layer로 정한다. 여러 절연층을 관통하는 via는 직접 지정한다.
  - `applicable_particle_source_steps`(CNO용, 13절)
- 입력 검증
  - contact 픽셀은 위·아래 layer 전극과 모두 겹쳐야 하며, 이어지는 두 전극은 같은 `net_id`여야 한다.
  - 한 contact hole이 같은 layer의 전극 둘에 걸치거나 어느 전극에도 닿지 않으면 오류다.
  - contact hole은 mask의 연결 성분 하나를 한 개로 본다.
- 자동 단자
  - 전극과 겹치는 contact hole 하나하나를 단자로 둔다.
  - ROI 경계에 닿은 전극 픽셀은 ROI 밖으로 이어지므로 단자로 둔다. 경계를 따라 이어진 한 구간을
    단자 하나로 본다.
  - 단자가 2개 미만인 전극은 판정하지 않고 `insufficient_terminals`로 기록한다.
- 판정(`src/particle_sim/failure/open_failure.py`)
  1. 입자가 없는 상태에서 단자를 연결 성분별로 묶는다(기준 분할). 이미 끊긴 단자는 처음부터 다른
     묶음이므로 입자 탓으로 세지 않는다.
  2. 이물 위치마다 기준 분할에서 같은 묶음이던 단자끼리 분리되면 open이다.
  3. ROI 경계 단자가 완전히 덮이면 open이다. contact 단자가 막히는 경우는 open에서 제외하고 CNO로
     계산한다. 막힌 contact 단자는 해당 위치의 open 연결 비교에서 빠진다.
  4. 단자가 없는 조각이 떨어지는 경우는 무시한다.
  5. 전극 bounding box 연결성 판정과 공정 간 결과 재사용을 그대로 적용한다.
- 적용 이물 유입 공정은 layer별 `open_applicable_steps_by_layer`로 명시한다. 적층 순서만으로
  자동 결정하지 않는다. 현재 파일럿 선언(M1의 G1이 `S_M1`, `S_M2`, `S_PLN`에 적용)이 "배선 형성 전
  낙하만 open"이라는 단순 규칙과 맞지 않아 물리 규칙을 코드에서 추정하면 안 되기 때문이다.
  schema 2.0에서 선언이 없으면 입력 오류다.
- 명시한 `open_paths`가 있는 전극은 그 경로가 자동 단자보다 우선한다. 기존 결과와 비교하기 위해서다.
  명시 경로에도 입자 없는 상태의 연결 검사를 적용해, 명목 입력에서 끊겨 있으면 입력 오류로 막고
  산포 조건에서 끊기면 `baseline_open`으로 따로 기록한다.
- 산포: 전극은 자기 layer의 overlay·CD로, contact는 `hole_layer`의 overlay·CD로 변형한다. 변형 후
  단자와 전극의 교집합이 비면 open으로 세지 않고 `terminal_lost`로 따로 기록한다.
- 출력: 불량 단위 ID `auto_open:<electrode_id>`, `unit_kind=open_electrode`, 판정에 쓴 단자 목록
  (contact ID, ROI 경계 구간)과 분류를 결과에 남긴다.

### 완료 기준

- 직선·우회·분기·ㄱ자 배선에서, 명시 단자 영역을 contact로 바꾼 자동 단자 결과가 기존 단자 방식과
  픽셀 단위로 일치한다(contact 막힘 위치는 CNO로 옮겨 간 것만 차이로 허용).
- 두 섬 전극, 산포로 끊긴 선, 단자 소실 사례에서 후보 전체가 open이 되지 않고 `baseline_open`,
  `terminal_lost`로 분류된다.
- 모서리 조각, 단자 없는 가지 절단은 open으로 세지 않는다.
- 전극 bounding box 판정이 전체 ROI 기준 구현과 모든 후보에서 일치한다.
- 200² 복합 파일럿의 TAT·RSS를 5회 재측정하고 기존 기록과 비교한다.

## 13. CNO (Contact Not Open) 불량 유형 추가

**2026-09-19 구현 완료(잠정 가림 기준 `full_cover_any_hole`).**

### 정의

CNO는 이물이 contact·via hole을 막아 위·아래 전극이 연결되지 않는 불량이다. 배선이 끊기는 open,
전극끼리 붙는 short와 다른 독립 불량 유형으로 계산한다. 보고서 표기는 `CNO`로 한다.

### 계산

- 불량 단위: contact 연결 하나(`cno:<contact_id>:<lower_electrode>__<upper_electrode>`),
  `failure_mode=contact_not_open`, `unit_kind=contact`
- hotspot: 이물 원판이 해당 연결의 contact hole을 막는 이물 중심 위치
  - 잠정 가림 기준: hole의 contact 픽셀이 **전부** 가려질 때 CNO. 일부 가림 기준은 11절 개선 항목이다.
  - 같은 연결에 hole이 여럿이면 하나라도 막히면 CNO(보수적 잠정 기준, 11절 개선 항목).
  - 계산은 hole별로 원판 가림 여부를 구하므로 연결성 labeling이 필요 없다. hole 주변 이물 반경
    범위만 후보로 계산한다.
- 적용 이물 유입 공정: contact별 `applicable_particle_source_steps`로 명시한다. schema 2.0에서 없으면
  입력 오류다. 추천 선언 범위는 아래 도체 형성 후부터 위 도체 증착 전까지 떨어진 공정이며,
  `TODO_PHYSICS`로 가정임을 기록한다.
- 산포: contact는 `hole_layer`의 overlay·CD로 변형한다. 변형 후 hole이 전극과 겹치지 않으면
  `contact_lost`로 따로 기록하고 CNO CA에 넣지 않는다.

### 출력과 위험도

- `ca.parquet`, `hotspots.zarr`, 확률 출력에 CNO 단위를 short·open과 같은 방식으로 추가한다.
- union은 `lateral_short`, `vertical_short`, `lateral_open`, `contact_not_open`별로 만들고, `union:any`는
  네 유형을 모두 OR한다. 한 위치에서 open과 CNO가 함께 생기면 `1/k` 배분으로 중복 없이 나눈다.
- 전극 기여도는 CNO를 위·아래 두 전극에 절반씩 배분한다(short와 동일). `risk_by_contact.parquet`를
  추가하고 `risk_summary.csv`에 `contact` 집계 수준을 넣는다.
- 기존 sim 비교에서 CNO와 자동 open 단위는 기존 계산 대상이 아니면 `new_unit`으로 분류한다(14절 5번과
  함께 수정).
- schema 2.0은 아직 확정 전이므로 새 열·값을 2.0 계약에 포함하고 문서의 열 설명을 같이 갱신한다.

### 완료 기준

- contact 완전 가림, 일부 가림, 가림 없음 합성 사례에서 잠정 기준대로 CNO가 판정된다.
- 적용 공정이 아닌 조건에서는 CNO CA가 0이고 `not_applicable_source_step`으로 분류된다.
- open·short·CNO가 겹치는 합성 사례에서 배분 CA 합과 `union:any` CA가 일치한다.
- Fab 위험도의 layer·전극·contact 집계 합이 union 기반 합과 일치한다.
- 결과 manifest와 보고서에 가림 기준과 중복 hole 기준이 잠정값임을 표시한다.

## 14. 2026-09-19 코드 검토 결함 (12·13절 외)

| 번호 | 우선순위 | 위치 | 결함 | 수정 방향 |
|---:|---|---|---|---|
| 3 | 완료 | `comparison/legacy_sim.py:88-130` | 비교 계약에 product·mask 식별자가 없고, run manifest에 pixel 크기가 없으면 검사를 건너뜀 | 계약에 product_id·design_version·mask hash를 추가해 `input_hashes`와 대조. pixel 값이 없으면 오류 |
| 4 | 완료 | `legacy_sim.py:181-198` | shape 검증 전에 출력 폴더와 zarr를 만들어 실패 시 폴더가 남고 재실행이 막힘 | 모든 reference 배열 검증 후 출력 생성 |
| 5 | 완료 | `legacy_sim.py:172-180, 233-256` | open 경로가 `new_pair`로 섞이고, 적용 공정이 아닌 open 조건에도 record를 요구 | `new_unit`·`not_applicable` 상태 분리 (단계 5a에서 반영) |
| 6 | 완료 | `pipeline.py:661, 745` | `convergence_delta`가 전체 seed·prefix의 최대값이 되고 `inf`를 버림 | 마지막 비교 기준으로 통일하거나 정의를 문서화하고 `inf` 보존 |
| 7 | 완료 | `aggregation/fab_weighting.py:131-148` | 평균 이물 표를 이전 `--ca` 경로에 넣으면 공정 ID와 Fab layer를 잘못 대조 | `mean_per_sample` 표는 `--run` 경로만 허용 |
| 8 | 완료 | `pipeline.py:131, 852` | 산포 임시 mask 폴더가 예외 시 출력 폴더에 남음 | `try/finally`로 정리 |

14절 결함은 2026-09-20에 모두 수정했다. 비교 계약은 schema 2.0으로 올리고 `product_id`,
`design_version`, `input_hashes`를 필수로 했다. `convergence_delta`는 seed별 마지막 prefix 비교의
최댓값이며 `inf`를 보존한다. 평균 이물 표는 `apply-fab --run`에서만 허용한다.

문서 수치 정정(반영 완료): `pair_margin_um` 크기는 파일 크기 합 621,304,169 B = 592.5 MiB다(594 MiB는 디스크
블록 사용량). `CODE-IMPLEMENTATION-PLAN.md`의 공정 간 cache 통합 RSS는 222MB다.

## 한 번에 수정할 때의 검증 묶음

1. 기존 단위·계약·CLI 테스트 전체를 실행한다.
2. 단순 합성 도형에서 기준 조건, 산포 0, 산포 64조건 결과를 비교한다.
3. 복합 200²와 1000²에서 전극 조합 누락, 크기 단조성, 중복 제거, CA 적분을 확인한다.
4. 겹치는 전극 조합에서 Fab 전극 조합별 단순합과 중복 제거 전체 가중의 차이를 확인한다.
5. `source_step != fab_distribution_layer`인 예제로 Fab와 simcompare 연결을 함께 검증한다.
6. open 최적화 전후를 모든 후보 픽셀에서 비교하고 TAT/RSS를 반복 측정한다.
7. 외부 기존 `sim` 동일 입력 비교에서 미분류 차이가 없는지 확인한다.
8. 자동 단자 open과 CNO의 합성 사례에서 기준 분할, 단자·contact 소실, 배분 합계를 확인한다.

최종 완료 조건은 코드 테스트 통과만이 아니다. 출력 schema, 입력·출력 가이드, simcompare 계약,
보고서 수치와 그림을 같은 버전으로 다시 생성해야 한다.
