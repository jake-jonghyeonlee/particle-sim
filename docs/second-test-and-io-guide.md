# 2nd 복합 형상 시험 계획과 입출력 가이드

이 문서는 `particle-sim`으로 기존 `sim`을 대체하기 전에 수행할 두 번째 로컬 시험의
범위와 제품 manifest YAML 작성법, 실행 결과의 의미를 설명한다. 설정 파일 안의 짧은 주석에
의존하지 않고도 입력을 작성하고 결과를 검토할 수 있게 하는 것이 목적이다.

실행 순서는 필수 코드 보완을 먼저 완료하고 수정 버전과 출력 schema를 고정한 뒤, 같은 입력으로
기존 `sim`과 hotspot·CA·TAT를 비교하는 것으로 확정한다. `sim-compare` 연결은 대체 판정 이후다.

YAML 문법 자체는 `#` 주석을 허용한다. 다만 주석은 파싱 결과와 `config_hash`에 포함되지 않고
실행 결과에도 전달되지 않는다. 따라서 물리값의 근거, 필드 의미, 승인 기준은 이 문서와 실행
manifest로 관리한다. 입력 YAML에는 계산에 필요한 값만 둔다.

보고서 본문과 설명 그림에서는 적층 순서를 쉽게 읽도록 `mask1`, `mask2`, `mask3`을 사용한다.
기존 합성 fixture의 실제 ID인 `M1`, `M2`, `M3`과 `S_M1` 등은 YAML·결과 파일에서만 유지한다.
새 실제 제품 manifest도 의미가 허용한다면 `mask1`, `mask2`처럼 바로 이해할 수 있는 ID를 권장한다.

검토 과정에서 확인된 코드 보완 항목과 수정 이유·완료 기준은
[`CODE-CHANGE-BACKLOG.md`](CODE-CHANGE-BACKLOG.md)에 별도로 정리했다.

## 1. 2nd test의 목표

첫 pilot은 200×200 px에서 4개 전극과 3개 전극 조합을 사용한다. M2에 A·B 두 전극이 있지만
형상이 대부분 직선이어서 실제 제품의 다중 전극과 다층 중첩을 충분히 대표하지 않는다.

두 번째 시험에서는 다음을 확인한다.

1. 한 레이어의 여러 전극을 서로 분리된 마스크로 정확히 읽는가.
2. 굴곡, 분기, 좁은 간격, 섬 형상을 가진 전극의 lateral short hotspot을 계산하는가.
3. 다른 레이어 전극이 정렬·부분 중첩·대각 교차할 때 vertical short hotspot 위치가 맞는가.
4. 여러 전극 조합의 hotspot이 겹칠 때 중복 제거 전체 CA가 픽셀 OR과 일치하는가.
5. 동일 입력으로 기존 `sim`과 hotspot·CA를 비교할 수 있는가.
6. 전극과 전극 조합 증가에 따라 기준 조건 계산 TAT, 최대 RSS, 출력 크기가 얼마나 증가하는가.

시험은 두 구간으로 분리한다. `2A`는 **기준 조건 short 검증**이며, open 계산은 후보 위치마다
연결성을 검사해 비용 특성이 다르므로 `2A` TAT에 섞지 않는다. `2B`는 같은 복합 fixture에
굴곡·분기 전극의 open 경로 3개를 추가해 lateral short, vertical short, open을 각각 단독 실행하고
통합 실행과 비교한다. `2A`와 `2B`의 200×200 합성 시험은 완료했으며, 산포 계산과 실제 제품 시험은
기준 조건 결과가 확인된 뒤 진행한다.

## 2. 시험 입력 설계

### 2.1 정확성 시험에서 pilot과 고정할 조건

입력 크기는 고정되어 있지 않다. 모든 마스크 shape가 같고 픽셀 크기로 계산한 실제 폭·높이가
ROI와 맞으면 1000×1000 이상도 사용할 수 있다. 200×200은 실제 운용 크기의 상한이 아니라,
형상 복잡도 외의 변수를 통제하고 픽셀 결과를 빠르게 검사하기 위한 정확성 fixture다.

정확성 시험에서는 아래 조건을 기존 pilot과 같게 유지한다.

| 항목 | 값 |
|---|---:|
| 정확성 fixture 이미지 | 200×200 px |
| 픽셀 크기 | x=0.5, y=0.5 µm/px |
| ROI | 100×100 µm |
| 도체 레이어 | M1, M2, M3 |
| 이물 유입 공정 | S_M1, S_M2, S_PLN |
| 이물 직경 | 2, 4, 6, 8, 10 µm |
| 매립 깊이/접촉 여유 | pilot과 동일한 0.0/0.0 µm |
| 후속 layer 정책 | pilot과 동일한 `protrusion`, 계수 1.5 |

정확성 검증을 통과한 뒤 같은 전극·전극 조합·이물 유입 공정·크기 구성을 1000×1000 px로 확장한 성능 fixture를
별도로 실행한다. 픽셀 크기를 0.5 µm/px로 유지하면 ROI는 500×500 µm가 된다. 실제 `sim`
입력이 더 크다면 1000×1000 결과를 확인한 다음 대표 실제품 크기를 세 번째 측정점으로 추가한다.

### 2.2 증가시킬 복잡도

fixture 이름은 `complex_multilayer`로 한다.

| 항목 | 기존 pilot | 2nd test |
|---|---:|---:|
| 전극 수 | 4 | 10 (M1=3, M2=4, M3=3) |
| 선언한 전극 조합 수 | 3 | 16 (lateral short=6, vertical short=10) |
| 이물 유입 공정 수 | 3 | 3 |
| 이물 크기 수 | 5 | 5 |
| open 경로 수 | 0 | 0 (`2B`에서 별도 추가) |

전극 형상에는 다음 요소를 각각 최소 한 번 포함한다.

- L자와 T자 분기
- 평행 구간과 국부적으로 좁아지는 간격
- 한 net에 속한 두 개의 떨어진 도형
- mask1↔mask2 및 mask2↔mask3의 완전 정렬, 부분 중첩, 대각 교차
- 의도적으로 멀리 떨어져 hotspot이 없어야 하는 음성 대조 전극 조합
- 세 전극 조합 이상의 hotspot이 겹쳐 동시 영향 조합 수가 1보다 큰 영역

같은 레이어의 전극 픽셀은 서로 겹치지 않아야 한다. 다른 레이어의 마스크는 겹쳐도 된다.
각 전극 마스크는 해당 `layer_mask` 안에 있어야 한다. 전극 외 TR 역할 영역이 있을 수 있으므로
전극 마스크 합집합이 `layer_mask` 전체와 같을 필요는 없다. 코드는 각 전극이 대응
`layer_mask`의 부분집합인지만 검사하고, 전극에 배정되지 않은 영역의 픽셀 수·면적·layer 대비 비율을
`validate-input`과 실행 manifest에 기록한다. TR은 전기 신호를 분리하는 역할이지만 현재 이물 해석에는
그 기능을 반영하지 않는다. 우선 남는 영역의 존재와 면적을 기록하고, TR의 신호 분리
효과를 불량 판정에 포함할지는 물리 기준과 검증 데이터가 준비된 뒤 결정한다.

### 2.3 전극 조합 선택 원칙

계산은 `electrode_pairs`에 선언한 조합만 수행한다. 다만 schema 2.0은 서로 다른 net의 모든 가능한
전극 조합을 감사하여, 계산 대상이면 `electrode_pairs`, 의도적 제외면 `excluded_electrode_pairs`에
이유와 함께 반드시 기록하게 한다. 어느 쪽에도 없는 조합이 하나라도 있으면 입력 검증이 실패한다.
`validate-input`은 선언·제외 상태와 함께 두 전극 mask의 최소 XY 거리, 적층 구간 사이의 최소 높이
차, 두 값을 결합한 최소 3D 거리를 출력한다. 이 값은 조합 검토용이며 이물 유입 공정·표면·후속 layer
정책을 적용한 hotspot 판정을 대신하지 않는다.
2nd test에서는 다음 규칙으로 16개를 고정하고 목록을 결과 보고서에 그대로 남긴다.

- 같은 레이어, 다른 net의 인접 조합 6개: `lateral_short`
- 다른 레이어, 다른 net이며 XY 도달 가능성이 있는 조합 8개: `vertical_short`
- hotspot이 없어야 하는 먼 거리 조합 2개: 음성 대조 `vertical_short`
- 같은 net 조합은 평가 대상에 넣지 않는다.

실제품 시험에서는 “눈에 보이는 인접 전극 조합”만 수동 선택하면 누락될 수 있다. 제품별 전극 조합 생성
규칙 또는 승인된 전극 조합 목록을 별도 입력으로 고정해야 한다.

## 3. 정답과 기존 sim 비교 방법

비교 전에 두 프로그램의 좌표계, 픽셀 크기, 직경/반경 정의, 경계 처리, 기존 방식 kernel 모양과
반복 반올림 규칙을 동일하게 맞춘다. 기존 `sim` 결과는 아래 공통 축으로 변환한다.

```text
source_step × particle_size_um × failure_unit_id × y × x
```

전극 조합·크기·이물 유입 공정별로 다음 값을 기록한다.

| 비교값 | 정의 |
|---|---|
| 기존 hotspot 픽셀 | 기존 `sim`의 binary 결과 |
| 내부 xy kernel 참고 픽셀 | `hotspots.zarr/legacy_failure`; 외부 기존 sim과 동일하다고 간주하지 않음 |
| 현재 2.5D 픽셀 | `hotspots.zarr/nominal_failure` |
| XOR 픽셀 | 두 binary map이 다른 위치 수 |
| 한쪽에만 있는 픽셀 | 기존 전용 / 현재 전용으로 분리 |
| CA 차이 | 절대 µm²와 기존 대비 상대 차이 |
| 위치 차이 | XOR overlay와 연결 성분별 bounding box |

`legacy_failure`는 xy kernel 차이를 설명하기 위한 내부 참고값이다. 외부 기존 sim도 높이를 반영한
2.5D 계산이므로 `legacy_failure`를 기존 sim 결과로 간주해서는 안 된다. 기존 sim과 현재 구현이
공통으로 계산하는 전극 조합은 좌표 반전, 경계 padding, kernel 모양과 거리 규칙을 맞춰 직접 비교한다.
G–A처럼 기존 sim의 계산 대상에 없던 대각 전극 조합은 기존값을 0으로 채우지 않고 `N/A` 또는
`new_pair`로 표시하며, 새 CA와 중복 제거 전체 CA 기여분을 별도로 제시한다.

fixture 자체의 독립적인 확인 조건도 둔다.

- 음성 대조 전극 조합의 hotspot은 0 px이다.
- 큰 직경 hotspot은 작은 직경 hotspot을 포함한다(고정된 기준 조건에서 확인).
- 각 전극 조합 결과는 두 전극 접촉 지도의 AND와 같다.
- `union_failure`는 소속 전극 조합 hotspot의 픽셀 OR와 같다.
- 중복 제거 전체 CA는 전극 조합별 CA 단순합보다 크지 않으며, 중첩 구간이 있으면 더 작다.
- 계산 제외가 승인된 전극 조합은 출력 축에 나타나지 않고 pair audit에는 `excluded`로 남는다.
- 다른 net인데 선언도 제외도 하지 않은 조합은 schema 2.0 입력 검증에서 거부된다.

## 4. TAT 측정 계획

### 4.1 측정 구간

TAT는 다음 세 구간을 분리한다.

1. 입력 검증: `validate-input`
2. 주 측정값: 기준 조건 계산 명령인 `simulate-nominal` 프로세스 시작부터 종료까지
3. 부가 측정값: `make-report`와 전체 `validate → simulate → report`

마스크 생성 시간과 사람이 전극 조합을 정의하는 시간은 계산 TAT에 포함하지 않는다. 기존 `sim`도
가능한 한 같은 범위로 계측한다.

### 4.2 실행 조건

- 같은 Mac, 같은 전원 모드, 다른 무거운 작업이 없는 상태에서 실행한다.
- 각 조합을 1회 워밍업한 뒤 5회 측정한다.
- 각 실행은 새 출력 디렉터리를 사용한다. 기존 결과 덮어쓰기는 허용되지 않는다.
- wall time, user/system CPU time, 최대 RSS, 결과 디렉터리 논리 크기를 JSON으로 보존한다.
- 대표값은 wall time 중앙값이며 최소/최대와 변동계수도 함께 표시한다.
- 입력 manifest와 모든 mask의 SHA-256, git commit, Python/패키지 버전을 기록한다.

### 4.3 비교표와 증가율

| 프로그램 | 입력 | 반복 | 측정 목적 |
|---|---|---:|---|
| 현재 particle-sim | pilot 200×200 | 5 | 현재 기준선 |
| 현재 particle-sim | complex 200×200 | 5 | 전극·전극 조합 복잡도 증가 비용 |
| 기존 sim | pilot/complex 200×200 | 각 5 | 픽셀 결과와 교체 TAT 기준 |
| 현재 particle-sim | complex 1000×1000 | 5 | 실제 운용 크기 성능 |
| 기존 sim | complex 1000×1000 | 5 | 동일 대형 입력의 교체 TAT 기준 |

보고할 핵심 값은 다음과 같다.

```text
복잡도 TAT 배수 = median(현재 complex) / median(현재 pilot)
교체 TAT 배수   = median(현재 complex) / median(기존 complex)
메모리 증가 배수 = peak_RSS(현재 complex) / peak_RSS(현재 pilot)
출력 증가 배수   = output_bytes(현재 complex) / output_bytes(현재 pilot)
```

동일 이미지 크기에서 전극 조합 조건 수는 pilot의 `이물 유입 공정 3 × 크기 5 × 전극 조합 3 = 45`에서
`3 × 5 × 16 = 240`으로 5.33배가 된다. 전극 contact 조건은 4개에서 10개로 2.5배다.
따라서 실행 시간은 대략 이 두 배수의 영향을 받지만, 이는 계획용 추정일 뿐 실제 결과로
대체해야 한다. 출력의 전극 조합별 지도 부분은 약 5.33배 증가할 것으로 예상된다. 200×200에서
1000×1000으로 늘리면 픽셀 수가 25배가 되므로, 이미지 크기와 전극 조합 증가를 합친 시험은 출력과
메모리 부담이 매우 커질 수 있다.

1000×1000, 이물 유입 공정 3개, 크기 5개, 전극 조합 16개 조건에서 Zarr 배열의 압축 전 논리 크기는 전극 조합별 결과와
중복 제거 유형 3개 기준으로 약 1.6 GB다. 실제 디스크 크기는 hotspot 압축률에 따라 달라진다.
디스크 여유 공간을 측정 전에 확인하고 압축 전 논리 크기와 실제 파일 크기를 함께 기록한다.

현재 저장된 참고 측정에는 1000×1000 px, 이물 유입 공정 3개, 전극 조합 3개, 크기 32개, open 제외 기준 조건 실행이
있다. 단일 실행에서 wall 16.07초, 최대 RSS 533.6 MB, 출력 파일 논리 크기 합계 3.35 MB였다.
이는 대형 입력이 실제 실행된 근거지만 반복 측정값이 아니며, 전극 조합 16개인 복합 fixture의 TAT를
대신하지 않는다(`docs/report/evidence/large-short.json`).

초기 시험에서는 TAT 합격 상한을 임의로 정하지 않는다. 측정 결과와 실제 제품의 ROI·전극·전극 조합
수를 바탕으로 운영 예산을 정한 뒤 성능 회귀 기준을 고정한다.

### 4.4 failure mode별 2B 파일럿 측정

같은 200×200 복합 형상에서 lateral short 6개 전극 조합, vertical short 10개 전극 조합,
open 3개 경로를 각각 단독 실행했다. open 경로는 G1, D1, D2의 서로 떨어진 두 단자 영역을 잇는
굴곡·분기 배선으로 정했다. 이물 유입 공정 3개와 이물 크기 5개를 사용했으며, open은 경로별로
지정한 적용 공정에서만 계산했다. 각 mode를 1회 워밍업한 뒤 순서를 회전하며 5회 측정했다.

| 실행 | 활성 계산 조건 | TAT 중앙값 | 단독 TAT 합 내 비율 | 조건당 시간 | 최대 RSS |
|---|---:|---:|---:|---:|---:|
| lateral short | 90 | 1.444초 | 7.31% | 0.0160초 | 198.0 MB |
| vertical short | 150 | 1.785초 | 9.04% | 0.0119초 | 205.7 MB |
| open | 35 | 16.523초 | 83.65% | 0.4721초 | 193.5 MB |
| 통합 | 275 | 17.862초 | 단독 합의 90.44% | 0.0650초 | 224.6 MB |

open의 총 TAT는 lateral short의 11.45배, vertical short의 9.26배였다. 활성 조건 하나당으로
정규화하면 각각 29.43배와 39.67배다. 원인은 현재 open이 이물 중심 후보마다 손상된 전체 전극의
연결 성분을 다시 계산하기 때문이다. 통합 실행은 공통 입력·접촉 계산을 공유하므로 세 단독 실행
TAT 합보다 약 9.6% 짧았다.

통합 결과에서 hotspot 중복 픽셀을 `1/k`로 나눈 CA 구성비는 lateral short 3.48%, vertical short
61.26%, open 35.27%였다. 이 비율의 합과 전체 중복 제거 CA는 모두 17,188 µm²로 일치했다.
이는 다섯 이물 크기와 적용 공정을 같은 비중으로 합한 **합성 형상의 기하 위험 구성비**다. Fab
크기별 평균 이물 개수로 가중한 실제 위험도나 불량률은 아니다. 이물 크기가 6, 8, 10 µm로
커질 때 open 비율은 각각 29.55%, 39.64%, 45.58%로 증가했다.

이 수치는 최적화 전 기준선이다. 원시 실행값은
`docs/report/evidence/benchmark-failure-modes-200.json`, 요약표는
`docs/report/evidence/benchmark-failure-modes-200.csv`에 있다.

### 4.5 동일 형상·크기 open 결과 재사용 후 측정

open의 연결성 결과는 같은 전극 형상, 단자, 이물 크기와 `contact_margin_um`이면 이물 유입 공정명에
따라 달라지지 않는다. 적용 가능한 공정인지 여부는 기존처럼 공정마다 판단하되, 명목 실행의 실제
연결성 계산은 open 경로·크기당 한 번만 수행하고 같은 형상의 다른 적용 공정에서 재사용하도록
수정했다. cache에는 binary open 결과만 저장한다. 산포 실행은 조건 수만큼 큰 지도가 누적되지 않도록
무제한 메모리 cache를 적용하지 않았다. 산포에서의 재사용은 bounded 또는 disk cache로 별도
구현해야 한다.

이 변경으로 2B fixture의 open 연결성 계산 횟수는 35회에서 15회로 줄었다. 동일 조건으로 1회
워밍업 후 5회 재측정한 결과는 다음과 같다.

| 실행 | 수정 전 TAT | 수정 후 TAT | 속도 향상 | 수정 후 최대 RSS |
|---|---:|---:|---:|---:|
| lateral short | 1.444초 | 1.436초 | 1.02배 | 198.3 MB |
| vertical short | 1.785초 | 1.760초 | 1.01배 | 207.4 MB |
| open | 16.523초 | 7.545초 | 2.19배 | 193.7 MB |
| 통합 | 17.862초 | 8.810초 | 2.03배 | 222.3 MB |

open TAT는 54.3%, 통합 TAT는 50.7% 감소했다. 수정 후에도 open은 단독 TAT 합의 70.25%를
차지하므로 다음 병목이다. 최적화 전후 `hotspots.zarr`의 5개 배열은 모든 픽셀이 같고 디렉터리
SHA-256도 `c9d7a12c...bb2ed`로 일치했다. `ca.parquet` 285행과 `ca_union.parquet` 60행은 실행마다
달라지는 `run_id`를 제외한 모든 값이 같았다. CA 구성비도 그대로 유지된다.

이 단계의 원시값은 `benchmark-failure-modes-200-after-open-cache.json`에 있다. 이후 전극 밖 빈
ROI를 연결 성분 판정에서 제외하는 두 번째 최적화를 이어서 수행했다.

### 4.6 전극 bounding box 연결성 판정 후 측정

기존 `open`은 후보마다 전체 200×200 ROI의 연결 성분을 계산했다. 전극 밖 픽셀은 항상 비도체이므로
연결 여부에 영향을 줄 수 없다. 이물 중심과 footprint는 원래 전체 좌표계를 유지하고, 손상 후 연결
성분 판정만 해당 전극의 최소 bounding box로 제한했다.

| open 경로 | 전체 ROI 픽셀 | 전극 bounding box 픽셀 | 전체 대비 |
|---|---:|---:|---:|
| G1_OPEN | 40,000 | 7,081 | 17.7% |
| D1_OPEN | 40,000 | 4,875 | 12.2% |
| D2_OPEN | 40,000 | 4,900 | 12.3% |

최종 5회 측정은 다음과 같다.

| 실행 | 최초 기준선 | 공정 간 cache | bounding box 적용 | 최초 대비 |
|---|---:|---:|---:|---:|
| lateral short | 1.444초 | 1.436초 | 1.447초 | 변화 없음 |
| vertical short | 1.785초 | 1.760초 | 1.761초 | 1.01배 |
| open | 16.523초 | 7.545초 | 3.135초 | **5.27배** |
| 통합 | 17.862초 | 8.810초 | 4.423초 | **4.04배** |

최초 기준선 대비 open TAT는 81.0%, 통합 TAT는 75.2% 감소했다. 최종 단독 TAT 합에서 mode별
비율은 lateral short 22.81%, vertical short 27.76%, open 49.43%다. 최종 최대 RSS는 open
193.7 MB, 통합 222.3 MB로 기준선과 유사하다.

최초 기준선과 최종 결과의 Zarr 5개 배열은 모든 픽셀이 같고 디렉터리 SHA-256도 일치했다. CA
285행과 union CA 60행도 `run_id`를 제외하면 완전히 같다. 작은 직선 배선에서는 전체 ROI 연결성
기준 계산과 bounding box 계산을 자동 테스트로 직접 대조한다. 최종 원시 측정은
`benchmark-failure-modes-200-after-open-crop.json`, 전체 동등성 근거는
`open-optimization-equivalence.json`에 있다.

이 결과는 국부 전극이 전체 ROI보다 작은 현재 합성 형상의 성능이다. 실제 전극이 ROI 대부분을
차지하거나 긴 우회 경로가 있으면 개선 폭이 작을 수 있다. 1000×1000은 동일 형상을 단순 타일링하면
한 전극 ID의 bounding box가 전체 ROI로 다시 커지므로, 실제 제품의 전극별 mask 구조를 대표하는
입력으로 측정해야 한다.

### 4.7 동일 형상 padded 1000×1000 확장성 측정

open 최적화 후 ROI 크기의 영향부터 분리하기 위해 200×200 복합 형상 자체는 바꾸지 않고,
1000×1000 canvas 중앙에 배치했다. 전극과 단자 좌표만 함께 이동했으므로 절대 CA와 union CA는
200×200 결과와 같다. 이 시험은 큰 빈 ROI와 결과 배열 저장 비용을 확인하는 시험이며, 긴 배선이나
실제 제품 layout의 대표 시험은 아니다.

| 실행 | padded 1000² TAT 중앙값 | 200² 대비 | 최대 RSS | 출력 크기 | 단독 TAT 합 내 비율 |
|---|---:|---:|---:|---:|---:|
| lateral short | 5.678초 | 3.92배 | 534.1 MB | 243.7 MB | 30.09% |
| vertical short | 8.522초 | 4.84배 | 650.3 MB | 377.8 MB | 45.16% |
| open | 4.669초 | 1.49배 | 333.4 MB | 0.07 MB | 24.74% |
| 통합 | 16.385초 | 3.70배 | 843.6 MB | 621.4 MB | — |

open은 후보와 연결 성분을 전극 bounding box 안에서 계산해 canvas의 빈 공간이 25배로 늘어도 TAT는
1.49배 증가했다. 반면 short의 `pair_margin_um`과 hotspot 결과는 전체 ROI 배열이므로 TAT·메모리·
출력이 더 크게 늘었다. 통합 출력 621.4MB 중 `pair_margin_um` 디렉터리가 621.3MB(594MiB)로
대부분을 차지했다. binary hotspot 배열들은 압축이 잘 되지만 float32 margin 지도는 빈 영역에서도
전체 전극 조합 축을 가진다.

Zarr 공간 chunk를 512에서 1024로 키운 단일 탐색 실행은 출력이 621.4MB에서 502.1MB로 줄었지만,
통합 TAT가 17.673초로 늘고 최대 RSS도 920.0MB로 증가했다. 따라서 기본 chunk는 512로 유지했다.
margin 지도의 선택 저장·정밀도·압축 정책은 출력 계약과 수치 의미를 먼저 정한 뒤 별도 최적화한다.

원시값은 `benchmark-failure-modes-padded-1000.json`, 200² 비교와 기각한 chunk 실험은
`padded-1000-analysis.json`, chunk 단일 실행은 `benchmark-zarr-chunk-1024-padded-1000.json`에
있다.

## 5. 시험 수행 순서와 완료 기준

### 단계 A — fixture와 입력 검증

1. `complex_multilayer` mask와 manifest를 생성한다.
2. 전극별 색상 overlay와 레이어별 합집합 그림을 만든다.
3. `validate-input`을 통과시키고 입력 hash를 보존한다.

완료 기준: 전극 겹침·누락 없이 전극 10개, 전극 조합 16개, 이물 유입 공정 3개, 크기 5개가 검증 결과에
표시되고 사람이 overlay를 확인한다.

### 단계 B — 현재 기준 조건 계산의 구조 검증

1. 현재 `particle-sim`의 기준 조건 계산을 실행한다.
2. 음성 대조, 크기 단조성, 전극 조합 AND, 중복 제거 OR를 자동 검사한다.
3. 대표 정렬·부분 중첩·대각 전극 조합의 hotspot 겹침 그림을 확인한다.

완료 기준: 위 불변조건이 모두 통과하고 설명할 수 없는 hotspot이 없다.

### 단계 C — 기존 sim 대조

1. 동일 mask와 파라미터로 기존 `sim`을 실행한다.
2. 결과를 공통 축과 binary dtype으로 변환한다.
3. 기존 방식 픽셀·CA 차이와 2.5D 변화 위치를 표와 이미지로 만든다.

완료 기준: 기존 방식 차이가 0이거나 모든 차이가 규칙 차이로 분류된다. 미분류 차이가 있으면
대체 판정을 중단하고 원인을 수정한다.

### 단계 D — TAT 측정

다음 순서로 크기를 올린다.

1. 200×200 pilot과 200×200 complex를 기존/현재에서 워밍업 1회 + 측정 5회 실행한다.
2. 정확성 검증을 통과하면 1000×1000 complex를 기존/현재에서 같은 방식으로 측정한다.
3. 실제 입력이 1000×1000보다 크면 대표 실제품 shape를 추가하되, 예상 출력 크기와 실행
   예산을 먼저 계산한다.

완료 기준: 원시 측정 JSON, 중앙값 비교표, TAT/RSS/출력 증가율, 실행 환경과 input hash가
남아 동일 장비에서 재실행할 수 있다.

### 단계 E — 후속 시험 결정

기준 조건 short, `2B` open 최적화와 padded 1000² 확장성 시험은 완료했다. 기존 `sim` 산출물이
확보되면 동일 입력 대조를 최우선으로 수행한다. 실제 제품형 1000² 이상 mask가 확보되면 긴 배선과
큰 전극 bounding box가 포함된 TAT를 별도로 측정하고, 이후 Sobol 64조건 산포 시험을 수행한다.
`sim-compare` 연결은 기존 `sim` 대체 결과와 출력 계약이 확정된 뒤 수행한다.

## 6. 입력 디렉터리 구성

manifest에 적는 파일 경로는 manifest 파일이 있는 디렉터리를 기준으로 한 상대 경로다.
절대 경로와 `..`로 바깥 디렉터리에 나가는 경로는 거부한다.

```text
case/
├── manifest.yaml
└── masks/
    ├── layer_M1.png
    ├── layer_M2.png
    ├── layer_M3.png
    ├── electrode_G1.png
    ├── electrode_G2.png
    ├── electrode_D1.png
    └── ...
```

마스크 입력 규칙:

- PNG 또는 TIFF의 2차원 이진 이미지여야 한다.
- 픽셀값은 0/1 또는 0/255만 허용한다.
- RGB, 회색조 중간값, resize·anti-alias된 이미지는 거부한다.
- 모든 layer/electrode mask의 shape는 같아야 한다.
- shape 자체에는 200×200이나 1000×1000 같은 고정 제한이 없다. 사용 가능한 메모리·디스크와
  실행 시간에 의해 실질적인 크기가 정해진다.
- `shape_x × pixel_size_x_um`과 ROI 폭, y 방향의 높이가 반 픽셀 허용오차 안에서 맞아야 한다.
- 전극 마스크는 비어 있을 수 없다.

## 7. 제품 manifest YAML 작성법

### 7.1 전체 구조 예시

아래 예시는 필드 배치를 설명한다. 실제 실행에서는 모든 patterned layer와 전극 파일이 실제로
존재해야 한다. `layer_M2.png`에는 D1·D2·D3 외 TR 역할 영역이 포함될 수 있으며, 각 전극
마스크는 `layer_M2.png`의 부분집합이어야 한다.

```yaml
schema_version: "2.0"
product_id: TEST_COMPLEX_02
product_family: TEST
design_version: V01
design_lineage_id: TEST_L1
fab_id: FAB_TEST
is_synthetic: true

coordinate_system:
  origin: upper_left
  x_direction: right
  y_direction: down
  pixel_size_x_um: 0.5
  pixel_size_y_um: 0.5
  roi_width_um: 100.0
  roi_height_um: 100.0

stack:
  - layer_name: M1
    process_order: 0
    material_type: conductor
    thickness_um: 0.3
    patterned: true
  - layer_name: GI
    process_order: 1
    material_type: inorganic
    thickness_um: 0.4
    planarization_factor: 0.0
    patterned: false
  - layer_name: M2
    process_order: 2
    material_type: conductor
    thickness_um: 0.3
    patterned: true
  - layer_name: PLN
    process_order: 3
    material_type: organic
    thickness_um: 2.0
    planarization_factor: 0.8
    patterned: false
  - layer_name: M3
    process_order: 4
    material_type: conductor
    thickness_um: 0.1
    patterned: true

layer_masks:
  - layer_name: M1
    path: masks/layer_M1.png
  - layer_name: M2
    path: masks/layer_M2.png
  - layer_name: M3
    path: masks/layer_M3.png

electrodes:
  - electrode_id: G1
    layer_name: M1
    net_id: NET_G1
    electrode_role: GATE
    path: masks/electrode_G1.png
  - electrode_id: D1
    layer_name: M2
    net_id: NET_D1
    electrode_role: DATA
    path: masks/electrode_D1.png
  - electrode_id: D2
    layer_name: M2
    net_id: NET_D2
    electrode_role: DATA
    path: masks/electrode_D2.png
  - electrode_id: D3
    layer_name: M2
    net_id: NET_D3
    electrode_role: DATA
    path: masks/electrode_D3.png
  - electrode_id: P1
    layer_name: M3
    net_id: NET_P1
    electrode_role: PIXEL
    path: masks/electrode_P1.png

particle_source_steps:
  - step_id: S_M1
    lands_after_layer: M1
    surface_model: conformal
    fab_distribution_layer: S_M1
  - step_id: S_M2
    lands_after_layer: M2
    surface_model: conformal
    fab_distribution_layer: S_M2
  - step_id: S_PLN
    lands_after_layer: PLN
    surface_model: flat
    fab_distribution_layer: S_PLN

electrode_pairs:
  - pair_id: D1__D2
    electrode_a: D1
    electrode_b: D2
    failure_mode: lateral_short
  - pair_id: D2__D3
    electrode_a: D2
    electrode_b: D3
    failure_mode: lateral_short
  - pair_id: G1__D1
    electrode_a: G1
    electrode_b: D1
    failure_mode: vertical_short
  - pair_id: D2__P1
    electrode_a: D2
    electrode_b: P1
    failure_mode: vertical_short

# 서로 다른 net이지만 이번 계산에서 제외한 모든 조합도 이유를 남긴다.
excluded_electrode_pairs:
  - {electrode_a: G1, electrode_b: D2, reason: "승인된 도달거리 밖"}
  - {electrode_a: G1, electrode_b: D3, reason: "승인된 도달거리 밖"}
  - {electrode_a: G1, electrode_b: P1, reason: "중간 절연 구조로 제외"}
  - {electrode_a: D1, electrode_b: D3, reason: "승인된 비인접 조합"}
  - {electrode_a: D1, electrode_b: P1, reason: "승인된 도달거리 밖"}
  - {electrode_a: D3, electrode_b: P1, reason: "승인된 도달거리 밖"}

open_paths: []

particle_sizes_um: [2.0, 4.0, 6.0, 8.0, 10.0]
embedding_depth_um: 0.0
contact_margin_um: 0.0

legacy_kernel:
  size: 5
  shape: square
  iterations: round

post_step_policy:
  mode: protrusion
  coverage_break_factor: 1.5
```

### 7.2 최상위 식별자

| 필드 | 필수 | 의미와 규칙 |
|---|---|---|
| `schema_version` | 예 | `"2.0"` 권장. `"1.0"`은 이전 입력 호환용 |
| `product_id` | 예 | 실행과 결과를 식별하는 제품 ID |
| `product_family` | 예 | Fab 분포 선택과 제품군 비교에 사용 |
| `design_version` | 예 | 설계 버전 |
| `design_lineage_id` | 예 | 설계 계보 ID |
| `fab_id` | 예 | 실제 Fab 분포를 적용할 때 선택하는 Fab ID |
| `is_synthetic` | 예 | 합성 fixture면 `true`, 실제 제품이면 `false` |

정의되지 않은 필드를 추가하면 오류가 발생한다. 이름 오타도 무시되지 않는다.

### 7.3 좌표계

| 필드 | 허용값/단위 | 설명 |
|---|---|---|
| `origin` | `upper_left` | 배열 `[0,0]`이 좌상단 |
| `x_direction` | `right` | 열 index가 증가하는 방향 |
| `y_direction` | `down` | 행 index가 증가하는 방향 |
| `pixel_size_x_um` | 양수 µm/px | x 방향 실제 크기 |
| `pixel_size_y_um` | 양수 µm/px | y 방향 실제 크기; x와 달라도 됨 |
| `roi_width_um` | 양수 µm | 이미지 폭과 일치해야 함 |
| `roi_height_um` | 양수 µm | 이미지 높이와 일치해야 함 |

코드의 배열 좌표는 `(y, x)`, 물리 좌표 표기는 `(x_um, y_um, z_um)`이다. 이전 open 단자 좌표를 사용할 때도 `(y, x)` 순서다. 실제 제품에는 해상도 변경에 안전한 단자 영역 mask를 권장한다.

### 7.4 stack

| 필드 | 필수 | 설명 |
|---|---|---|
| `layer_name` | 예 | stack 안에서 유일한 레이어 이름 |
| `process_order` | 예 | 0 이상의 고유한 형성 순서 |
| `material_type` | 예 | `conductor`, `inorganic`, `organic`, `other` |
| `thickness_um` | 예 | 양수 두께. z는 process order대로 누적 |
| `patterned` | 예 | `true`면 `layer_masks`에 정확히 하나 필요 |
| `planarization_factor` | 조건부 | 0~1. `partial` 표면 또는 `protrusion` 정책의 비도체층에 필요 |
| `z_bottom_um` | 아니오 | 자동 누적 z와의 교차검사용; 값을 주면 일치해야 함 |
| `thickness_sigma_um` | 아니오 | 현재 명목 계산에서는 사용하지 않음 |

한 stack layer에 속한 모든 전극은 같은 명목 두께와 z 구간을 공유한다. 같은 이름의 레이어
안에 높이·재질·형성 순서가 다른 도체가 있다면 stack layer를 나눠야 한다.

### 7.5 layer_masks와 electrodes

`layer_masks`는 patterned layer의 전체 패턴이다. 레이어당 정확히 하나다. `electrodes`는 그
레이어 패턴 중 현재 short 평가 대상으로 구분한 전극 목록이므로 한 레이어에 여러 개를 둘 수 있다.
TR 영역은 이 목록에 포함하지 않으며 현재 불량 판정에도 사용하지 않는다.

| 관계 | 검증 규칙 |
|---|---|
| 전극 → stack | 전극 layer는 `patterned: true`인 conductor여야 함 |
| 같은 layer의 전극 | 픽셀 중첩 금지 |
| 전극 → layer mask | 각 전극 mask가 대응하는 layer mask의 부분집합이어야 함 |
| layer mask의 남는 영역 | TR 또는 비평가 영역으로 허용하고 픽셀 수·면적 비율을 기록. TR의 전기적 기능은 현재 미반영 |
| 전극이 없는 patterned layer | topography 전용 layer로 허용 |
| 다른 layer의 전극 | XY 픽셀 중첩 허용 |
| `net_id` | 전기적으로 같은 net을 나타냄; 같은 net끼리는 short 전극 조합으로 지정할 수 없음 |

하나의 `electrode_id` 마스크에 떨어진 여러 영역을 넣을 수 있다. 이 경우 short 계산은 전체를
하나의 전극으로 취급한다. open 연결성을 계산할 때는 단자와 ROI 내 연결 상태가 의미에 맞는지
별도로 확인해야 한다.

### 7.6 particle_source_steps

이 항목은 이물 유입 공정을 정의한다. YAML 필드명은 코드 계약이므로 그대로 사용한다.

| 필드 | 설명 |
|---|---|
| `step_id` | 결과의 이물 유입 공정 축에 기록되는 고유 ID |
| `lands_after_layer` | 이물이 유입되는 시점 직전에 형성된 stack layer |
| `surface_model` | `flat`, `conformal`, `partial` 중 하나 |
| `fab_distribution_layer` | 후속 Fab 이물 표의 측정 공정과 연결할 이름 |

- `flat`: landing layer의 평균 상면 높이를 전체 ROI에 사용한다.
- `conformal`: 아래 patterned layer의 요철을 그대로 반영한다.
- `partial`: 위 비도체층의 `planarization_factor`로 남는 요철을 계산한다.

### 7.7 electrode_pairs

이 항목은 평가할 전극 조합을 정의한다. YAML 필드명은 코드 계약이므로 그대로 사용한다.

| 필드 | 설명 |
|---|---|
| `pair_id` | 전극 조합을 식별하며 결과 축과 CA 표에 기록되는 고유 ID |
| `electrode_a`, `electrode_b` | `electrodes`에 존재하는 서로 다른 ID |
| `failure_mode` | `lateral_short` 또는 `vertical_short` |

`lateral_short`는 같은 layer의 lateral short, `vertical_short`는 다른 layer의 vertical short를
뜻하며 두 전극의 `net_id`는 달라야 한다. schema 2.0에서는 다른 net의 모든 조합이
`electrode_pairs` 또는 `excluded_electrode_pairs` 중 정확히 한 곳에 있어야 한다. 제외 항목은
`electrode_a`, `electrode_b`, 구체적인 `reason`을 갖는다.

### 7.8 open_paths

open은 이물이 배선을 가린 뒤에도 두 단자 영역 사이에 도전 경로가 남는지 판정한다. 실제 제품에는
두 단자 영역 mask와 open을 계산할 이물 유입 공정을 명시한다.

```yaml
open_paths:
  - path_id: data1_open
    electrode: D1
    terminal_a_mask_path: masks/D1_terminal_a.png
    terminal_b_mask_path: masks/D1_terminal_b.png
    applicable_particle_source_steps: [S_M2, S_PLN]
    fallback: error
```

두 단자 mask는 전체 입력과 같은 shape의 이진 이미지여야 하며, 비어 있지 않고 대상 전극의
부분집합이고 서로 겹치지 않아야 한다. 이물이 단자 영역 전체를 덮거나 손상된 전극에서 두 단자 영역을
잇는 4-연결 경로가 사라지면 open이다. 목록에 없는 이물 유입 공정에서는 해당 open 결과를 0으로
기록해 출력 축을 안정적으로 유지한다.

`terminal_a_yx_px`와 `terminal_b_yx_px` 좌표 쌍은 이전 입력 호환용이며 단자 mask와 함께 쓸 수 없다.
단자가 없을 때 `fallback: roi_boundary`를 쓰면 결과가 `open_proxy`로 표시된다. schema 2.0에서는
`applicable_particle_source_steps`가 비어 있으면 거부한다. `skeleton_cut`은 지원하지 않는다.

### 7.8.1 contacts와 자동 단자 open, CNO

contact·via 원본 mask를 넣으면 전극마다 open 단자를 자동으로 만들고, contact가 막히는 CNO(Contact Not
Open)를 별도 불량 유형으로 계산한다. schema 2.0에서만 쓴다.

```yaml
contacts:
  - contact_id: CNT_GI
    path: masks/contact_GI.png
    hole_layer: GI              # 구멍이 뚫리는 절연층
    # lower_layer / upper_layer: 생략하면 stack process_order에서 바로 아래·위 patterned 도체
    applicable_particle_source_steps: [S_M1]   # CNO를 계산할 이물 유입 공정
open_applicable_steps_by_layer:
  M1: []                        # contact로 이어졌지만 자동 open은 계산하지 않음
  M2: [S_M2, S_PLN]             # M2 전극마다 자동 단자 open 계산
```

- contact mask의 연결 성분 하나가 hole 하나다. hole은 아래·위 layer에서 각각 정확히 한 전극과 겹쳐야
  하고 두 전극은 같은 `net_id`여야 한다. 실제 접촉 픽셀은 hole ∩ 아래 전극 ∩ 위 전극이다.
- contact로 이어진 layer는 모두 `open_applicable_steps_by_layer`에 적는다. 빈 목록은 그 layer의 자동
  open을 계산하지 않는다는 뜻이다. 적용 공정은 적층 순서로 추정하지 않고 명시한다.
- **자동 단자 open(`auto_open:<전극>`):** 전극에 닿은 contact hole과 ROI 경계에 닿은 구간이 단자다.
  이물 없는 상태에서 이어져 있던 단자끼리 분리되거나 ROI 경계 구간이 완전히 덮이면 open이다. 단자가
  없는 조각이 떨어지는 경우는 무시한다. 단자가 2개 미만이면 `insufficient_terminals`로 기록한다.
  같은 전극에 `open_paths`가 있으면 그 경로가 우선한다.
- **CNO(`cno:<contact>:<아래 전극>__<위 전극>`):** 잠정 기준은 이물이 hole의 실제 접촉 픽셀을 모두
  덮으면 CNO, 같은 두 전극을 잇는 hole 중 하나라도 막히면 CNO다. 일부 가림과 병렬 hole 기준은
  논의 중인 개선 항목이며 실행 manifest의 `cno_coverage_rule`과 경고에 표시된다.
- contact가 막히는 위치는 open이 아니라 CNO로 계산한다. 두 유형이 같은 위치에 겹치면 union과 `1/k`
  배분이 중복 없이 처리한다.
- 산포에서 contact는 `hole_layer`의 overlay·CD를 따른다. 산포 YAML의 overlay·CD 키에 hole layer 이름을
  쓸 수 있다. 변형 후 전극과 이어지지 않는 hole은 `contact_lost`로 기록한다.
- open 경로는 이물 없는 상태에서 두 단자가 이어져 있어야 한다. 명목 입력에서 끊겨 있으면 입력 오류이고,
  산포 조건에서 끊기거나 단자가 사라지면 `baseline_open`·`terminal_lost`로 기록해 CA에서 제외한다.

### 7.9 이물과 후속 layer 정책

| 필드 | 설명 |
|---|---|
| `particle_sizes_um` | 양수이며 엄격히 증가하는 직경 목록 |
| `embedding_depth_um` | 0 이상. 이물 중심 높이 계산에 사용 |
| `contact_margin_um` | 현재 시험은 0 µm로 고정. 비영값은 이물 재질·계면의 실측 또는 물리 근거 확보 후 사용 |
| `legacy_kernel` | 기존 방식의 홀수 kernel 크기, `square/disk`, 반복 반올림 규칙 |
| `post_step_policy.mode` | `exclude`, `protrusion`, `wrap_contact` |
| `coverage_break_factor` | `protrusion`에서 필수인 양수 피복 파괴 계수 |

현재 `particle_sizes_um`은 원형 footprint와 구형 높이 모델의 직경으로 해석한다. 타원형·길쭉한
형상·불규칙 형상과 방향각에 따른 결과 차이는 현재 범위에 포함하지 않으며 추후 형상 모델 검토
항목으로 남긴다.

- `exclude`: 이물 유입 공정 이후 형성되는 전극과의 접촉을 제외한다.
- `protrusion`: 돌출 높이와 평탄화·피복 가정으로 접촉 가능성을 정한다.
- `wrap_contact`: 이후 전극의 접촉을 보수적으로 허용한다.

이 값들은 실제 제품 근거 없이 기본값처럼 재사용하면 안 된다. 실행 결과의
`run_manifest.json`에 `TODO_PHYSICS` 경고로 남는다.

### 7.10 산포 YAML과 Fab 실측 overlay·CD

제품 manifest는 기준 형상과 적층을 정의하고 산포 YAML은 mask layer별 Fab 실측 overlay·CD를
정의한다. schema 2.0에서는 x/y를 따로 쓰고 기준 layer, 부호, 자료 출처를 명시한다.

```yaml
schema_version: "2.0"
assumption_set: fab_overlay_cd_2026Q3
evidence: measured
overlay_reference_layer: M1
overlay_sign_convention: positive_x_right_positive_y_down
independence_model: independent
overlay_um:
  M2:
    x: {mean: 0.0, sigma: 0.15, evidence: measured, source_id: OVL-2026Q3-M2-X}
    y: {mean: 0.0, sigma: 0.18, evidence: measured, source_id: OVL-2026Q3-M2-Y}
cd_bias_um:
  M2:
    mean: 0.0
    sigma: 0.05
    evidence: measured
    source_id: CD-2026Q3-M2
```

숫자는 형식 설명용이다. `source_hash`, `measurement_period`, `fab_process`, `product_family`도 넣을 수
있다. `evidence: measured`인 overlay/CD 항목은 `source_id`가 필수다. 현재 계약은 mask별 CD와
layer별 overlay, x/y를 모두 독립으로 표본화하며 `independence_model`도 `independent`만 허용한다.
상관관계를 쓰게 되면 별도 schema가 필요하다.

overlay는 기준 layer 대비 상대 이동이다. 양의 x는 오른쪽, 양의 y는 아래쪽이다. 같은 layer의
`layer_mask`, 모든 전극 mask, 단자 mask에는 같은 overlay가 적용되고 단자 mask 자체에는 CD 확장을
적용하지 않는다.

`cd_bias_um`은 한쪽 mask 경계의 이동 거리다. 양수는 확장, 음수는 축소다. Fab 원자료가 전체 선폭
변화 `ΔCD`라면 대칭 변화 가정에서 입력은 보통 `ΔCD/2`다. 원자료가 이미 edge 이동량이면 다시
나누지 않는다. 현재는 정규분포만 지원하며 원본 표본·비정규분포는 후속 범위다.

### 7.11 Fab layer별 평균 이물 개수 표

보유한 Fab 자료는 layer별·이물 크기 구간별 sample당 평균 이물 개수이며, 이물 시뮬레이션 완료 후
적용한다. CSV 또는 Parquet에 다음 열을 둔다.

| 필드 | 필수 | 의미 |
|---|---|---|
| `fab_id`, `product_family` | 예 | 적용할 Fab와 제품군. 공통 분포는 제품군 `*` 사용 가능 |
| `fab_distribution_layer` | 예 | manifest의 동일 필드와 연결되는 Fab 측정 layer |
| `size_bin_left_um`, `size_bin_right_um` | 예 | 이물 크기 구간 |
| `size_representative_um` | 예 | `exact_point`로 연결할 시뮬레이션 직경 |
| `mean_particle_count` | 예 | 검사 sample 하나당 평균 이물 개수 |
| `inspection_area_per_sample_mm2` | 예 | sample 하나의 검사 면적 |
| `mean_density_per_mm2` | 아니오 | 비우면 평균 개수/검사 면적으로 계산; 값이 있으면 식과 일치해야 함 |
| `sample_count` | 아니오 | 평균에 사용한 sample 수 |
| `measurement_period`, `particle_type` | 예 | 측정 기간과 이물 종류 |
| `source_id`, `source_hash` | 아니오 | Fab 원자료 식별자와 hash |

```text
mean_density_per_mm2 = mean_particle_count / inspection_area_per_sample_mm2
```

평균 개수 입력은 현재 `--bin-method exact_point`만 허용한다. `fab_distribution_layer`에 적용하는 전체
위험은 그 공정에서 선언된 모든 `lateral_short`, `vertical_short`, `open` hotspot을 OR한 union CA에
밀도를 곱한 값이다. 이는 해당 layer 전극 mask의 합집합 면적이 아니다.

## 8. 실행 순서

저장소 루트에서 실행한다.

```bash
uv run particle-sim validate-input \
  --product-config work/complex_multilayer/manifest.yaml

uv run particle-sim simulate-nominal \
  --product-config work/complex_multilayer/manifest.yaml \
  --output work/complex_multilayer-run

uv run particle-sim make-report \
  --run work/complex_multilayer-run \
  --product-config work/complex_multilayer/manifest.yaml
```

출력 디렉터리는 존재하지 않거나 비어 있어야 한다. 기존 결과를 덮어쓰지 않는다.

산포 시험은 기준 조건 검증 뒤 별도 디렉터리에서 실행한다.

```bash
uv run particle-sim simulate-probability \
  --product-config work/complex_multilayer/manifest.yaml \
  --scenario-config configs/scenarios/nominal.yaml \
  --sobol-m 6 \
  --seed 42 \
  --additional-seed 43 \
  --additional-seed 44 \
  --output work/complex_multilayer-prob
```

`sobol-m 6`은 산포가 0이 아닐 때 seed당 64조건을 뜻한다. 위 예시는 세 seed, 총 192조건이다.
시나리오 YAML의 `evidence`가 `assumed`이면 결과는 가정한 산포에 조건부이며 실제 불량확률이 아니다.

기존 sim 결과와 기준 조건을 비교할 때는 `configs/legacy-sim-reference.example.yaml`을 복사해 실제
파일 경로와 모든 공통 조건을 채운 뒤 실행한다.

```bash
uv run particle-sim compare-sim \
  --reference-config work/legacy-reference.yaml \
  --run work/complex_multilayer-run \
  --output work/complex-legacy-comparison
```

Fab 평균 이물 개수를 적용해 layer·전극 조합·open·개별 전극 관점을 함께 생성한다.

```bash
uv run particle-sim apply-fab \
  --run work/complex_multilayer-prob \
  --fab-data work/fab-mean-particles.parquet \
  --product-config work/complex_multilayer/manifest.yaml \
  --bin-method exact_point --mode actual \
  --output work/complex-risk
```

## 9. 기준 조건 출력 설명

### 9.0 Zarr와 Parquet를 왜 나누어 쓰는가

두 형식은 같은 결과를 서로 다른 방법으로 보여준다.

- **Zarr는 지도 묶음이다.** 조건마다 생성되는 2차원 hotspot 지도를 여러 장 쌓아 저장한다. 예를 들어
  `이물 유입 공정 × 이물 크기 × 전극 조합 × y × x`처럼 5차원 배열을 보존한다. 큰 지도를 작은
  chunk로 나누므로 1000×1000 이상의 결과도 필요한 조건과 위치만 읽을 수 있다. 일반 문서 파일이나
  Excel 표는 아니며 Python으로 읽거나 PNG 미리보기로 확인한다.
- **Parquet는 분석용 표다.** Excel처럼 행과 열이 있지만, 필요한 열만 빠르게 읽고 압축하기 좋은
  이진 형식이다. 한 행에 Fab layer, 이물 크기, 전극 조합, CA와 단위를 저장하므로 정렬·필터·집계에
  적합하다. 사람이 바로 열어 볼 핵심 순위는 별도 CSV로도 내보낸다.

Zarr에서 “어느 위치가 위험한가”를 확인하고, Parquet에서 “어느 조건·전극의 숫자가 큰가”를 확인한다.
Parquet 값이 이상하면 같은 조건의 Zarr 지도로 돌아가 위치와 중첩을 점검한다.

### 9.1 hotspots.zarr

모든 배열의 공간 축은 `(y, x)`다.

| 배열 | dtype과 축 | 의미 |
|---|---|---|
| `legacy_failure` | uint8 `[step,size,unit,y,x]` | 기존 Erode/Dilate 방식의 0/1 hotspot |
| `nominal_failure` | uint8 `[step,size,unit,y,x]` | 기준 조건 2.5D의 0/1 hotspot |
| `pair_margin_um` | float32 `[step,size,unit,y,x]` | 두 전극 중 더 작은 접촉 여유 |
| `union_failure` | uint8 `[step,size,mode,y,x]` | 불량 유형별 및 전체 hotspot의 픽셀 OR |
| `union_pair_count` | uint16 `[step,size,mode,y,x]` | 해당 위치에서 동시에 영향받은 전극 조합 수 |

실제 축 이름은 Zarr 속성 `source_steps`, `particle_sizes_um`, `failure_unit_ids`,
`unit_kinds`, `union_keys`에서 읽는다. `pair_ids`는 이전 reader 호환용 별칭이다. 숫자 위치를 문서 밖에서 고정하지 않는다.

### 9.2 ca.parquet

한 행은 이물 유입 공정 × 이물 크기 × 불량 단위다. 실제 식별자는 `failure_unit_id`이며
`unit_kind`가 `electrode_pair`이면 `electrode_pair_id`, `open_path`·`open_electrode`이면
`open_path_id`, `contact`(CNO)이면 `contact_id`가 채워진다. CNO의 `electrode_a`·`electrode_b`는 아래·위
전극이다.

- `ca_legacy_um2`: 기존 방식 hotspot 면적
- `ca_nominal_um2`, `ca_nominal_mm2`: 기준 조건 2.5D 면적
- `ca_nominal_normalized`: ROI 면적 대비 비율
- `classification`: 직접 기하 접촉, potential vertical short, open 구분
- `ca_mean_um2`, `ca_p05_um2`, `ca_p50_um2`, `ca_p95_um2`: 산포 조건별 CA 통계
- `ca_allocated_*`: 겹친 픽셀을 활성 불량 단위 수 `k`로 나눈 `1/k` 배분 CA
- `ca_unique_*`: 해당 불량 단위만 담당하는 고유 픽셀 CA

### 9.3 ca_union.parquet

같은 조건의 전극 조합별 hotspot을 OR해 겹친 위치를 한 번만 센 결과다.

- `ca_union_nominal_um2`: 기준 조건의 중복 제거 전체 CA
- `ca_pair_sum_nominal_um2`: 전극 조합별 CA의 단순합
- 두 값의 차이: 여러 전극 조합의 hotspot이 같은 픽셀에 중첩된 면적 효과
- `failure_mode=any`: lateral short·vertical short·open·CNO를 모두 합친 중복 제거 결과

이 union은 전극 mask 자체를 합친 면적이 아니다. 한 `source_step`과 이물 크기에서 선언된 불량
hotspot의 union이다. 예를 들어 mask1 형성 후 들어온 이물이 mask1 내부 lateral short와
mask1–mask2 vertical short를 만들 수 있으면 두 결과가 함께 포함된다. Fab layer별 평균 이물 개수는
그 layer에 대응하는 `source_step`의 `failure_mode=any` union CA에 적용한다.

원인 분석용으로 lateral short의 전극 layer별 union과 vertical short의 layer 조합별 union을 추가
출력할 수 있지만, 이들을 다시 단순 합산해 전체값으로 사용하지 않는다.

### 9.4 기타 파일

| 파일 | 내용 |
|---|---|
| `legacy_vs_edt.csv` | 내부 기존 방식과 EDT 기준 결과의 픽셀 차이 및 원인 분류 |
| `run_manifest.json` | 입력 hash, 설정 hash, 단위, 버전, 물리 가정 경고, 실행 ID |
| `report/report.md` 및 PNG | 사람이 검토할 layout, 전극 조합별 hotspot, CA, 중복 제거 그림 |

### 9.5 Fab 적용 위험도 출력

최종 위험도는 하나의 합계로 축약하지 않고 다음 세 관점으로 출력한다.

| 관점 | 기본 분석 단위 | 대표 출력 | 해석 |
|---|---|---|---|
| Fab layer별 | `fab_distribution_layer × particle_size` | `risk_by_layer.parquet` | 해당 공정에서 유입된 이물의 전체 중복 제거 위험 |
| 전극 조합별 | `electrode_pair_id × particle_size` | `risk_by_pair.parquet` | 어느 실제 전극 조합이 민감하고 전체에 얼마나 기여하는지 |
| layer와 전극 동시 | `fab_distribution_layer × electrode_pair_id × particle_size` | `risk_by_pair.parquet` | 어느 공정 이물이 어느 전극 조합에 크게 연결되는지 |
| 개별 전극별 | `electrode_id` | `risk_by_electrode.parquet` | 여러 조합에 반복 등장하는 주요 전극 |

`risk_by_pair.parquet`에는 `electrode_pair_id`, `electrode_a`, `electrode_b`, 각 전극 layer,
`failure_mode`, Fab layer와 다음 지표를 함께 둔다.

- `raw_pair_risk`: 해당 전극 조합의 전체 hotspot CA에 Fab 밀도를 적용한 값이다. 전극 조합 자체의
  민감도를 비교하는 데 사용하며 여러 조합의 값을 더하지 않는다.
- `allocated_pair_contribution`: 같은 픽셀·산포 조건에서 `k`개 조합이 동시에 hotspot이면 각 조합에
  `1/k`를 배분한 뒤 Fab 밀도를 적용한 값이다. 같은 Fab layer와 크기에서 전극 조합별 값을 합하면
  `risk_by_layer.parquet`의 union 위험과 일치한다.
- `unique_pair_contribution`: 해당 전극 조합을 제거했을 때 union이 줄어드는 양이다. 다른 조합과 완전히
  겹치면 0이 될 수 있으므로 이 값만으로 중요도를 판단하지 않는다.

산포 결과는 평균 확률 지도에서 나중에 나누지 않는다. 각 산포 조건에서 먼저 `1/k`를 배분하고 그
결과를 조건 전체에서 평균해야 전극 조합별 합계와 union 평균이 일치한다.

개별 전극 출력은 두 해석을 구분한다.

- `involvement_risk`: 그 전극이 포함된 전극 조합의 배분 기여를 모두 모은 값이다. “자주 관여하는
  전극”을 찾는 순위이며 전극 전체의 합은 union과 일치하지 않을 수 있다.
- `allocated_electrode_contribution`: short 조합의 기여를 두 전극에 절반씩 나누고 open은 해당 전극에
  전부 배분한 값이다. 개별 전극 전체의 합이 union 위험과 일치한다.

`1/k`와 절반 배분은 전기적 고장 원인을 판정하는 물리식이 아니다. 겹치는 hotspot을 중복 없이 설명하기
위한 중립적인 보고 규칙이며, 원본 `raw_pair_risk`도 항상 함께 보존한다. 사람이 바로 확인할 수 있도록
각 관점의 상위 항목을 `risk_summary.csv`에 출력하고, 단위·배분 규칙·Fab 입력 의미·출력 hash는
`risk_manifest.json`에 보존한다. open은 `risk_by_open_path.parquet`으로 따로 확인한다.
위험도 실행은 원본 `ca.parquet`·`ca_union.parquet`의 hash와 run/product ID를 먼저 검증하며, 제품
YAML·Fab 표·원본 run manifest의 SHA-256도 `risk_manifest.json`에 남긴다.

## 10. 산포 출력 설명

| 파일/배열 | 의미 |
|---|---|
| `probability.zarr/probability` | 조건별 0/1 hotspot을 평균한 `[0,1]` 지도 |
| `contact_margin_mean_um` | 조건별 접촉 여유 평균 |
| `union_probability` | 각 조건에서 먼저 OR한 뒤 평균한 중복 제거 전체 확률 지도 |
| `ca.parquet` | 기준 조건 CA, 산포 평균 CA, 조건별 CA p05/p50/p95, 조건 수 |
| `ca_union.parquet` | 불량 유형별 중복 제거 확률 CA; 실행 종류와 함께 해석 |
| `convergence.csv` | 64/128/256조건 누적 구간의 CA와 직전 대비 변화 |
| `run_manifest.json` | seed, Sobol m, scenario 근거와 경고 |

확률 지도는 입력한 공정 산포 모델에 조건부인 **기하 접촉 위험**이다. 전기적 불량확률이나
Fab 수율로 직접 해석하지 않는다.

Sobol은 먼저 0~1 누적확률 공간을 고르게 채우고, 각 좌표를 지정한 분포의 역누적분포로 변환한다.
따라서 정규분포를 지정하면 최종 변화량은 정규분포를 근사한다. 다만 유한 조건의 표본 통계는 입력
평균·표준편차와 정확히 같지 않을 수 있고, 실제 Fab 분포가 비정규이거나 변수 간 상관이 있으면 해당
분포와 상관구조를 별도로 입력해야 한다.

수렴은 동일 seed의 중첩 표본 64→128→256으로 확인한다. 중복 제거 전체 평균 CA는 128→256 상대
변화 1% 이하 또는 절대 변화 1픽셀 면적 이하, 주요 전극 조합은 2% 이하 또는 1픽셀 면적 이하를
초기 기준으로 권장한다. 기준을 넘으면 512조건으로 늘린다. 최종 pilot은 선택한 조건 수에서 Sobol
scramble seed를 3개 이상 바꾸어 중복 제거 전체 평균 CA의 최댓값과 최솟값 차이가 seed 평균의
1% 이내인지 추가 확인한다. 평균 CA가 0에 가까우면 1픽셀 면적의 절대 기준을 사용한다.

Fab에서 작은 이물이 많고 큰 이물이 적은 관측은 이물 크기별 발생 분포로 처리한다. 각 크기의 CA를
먼저 계산한 뒤 크기 구간별 이물 밀도로 가중하며, 이 분포를 Sobol 공정 산포 변수로 섞지 않는다.

## 11. 결과 검토 체크리스트

- [x] `validate-input` 출력의 이미지 크기, 전극 수, 전극 조합 수, 이물 유입 공정이 계획과 같다.
- [x] 레이어 overlay에서 전극 누락·동일 레이어 중첩·축 반전이 없다.
- [x] Zarr attribute의 축 이름과 manifest 순서가 맞는다.
- [x] 대표 전극 조합의 0/1 hotspot과 CA를 픽셀 면적으로 역산할 수 있다.
- [x] 중복 제거 전체 지도가 전극 조합별 지도의 OR와 같고, 전체 CA가 단순합보다 크지 않다.
- [ ] 기존 `sim` 비교에서 미분류 XOR 픽셀이 없다.
- [x] TAT 결과에 중앙값, 변동폭, 최대 RSS, 출력 크기, 입력 hash가 있다.
- [x] `run_manifest.json`의 `TODO_PHYSICS` 값을 실제 시험 가정과 대조했다.

## 12. 2026-09-18 수행 결과

복합 합성 fixture의 현재 `particle-sim` 기준 조건 시험을 완료했다. fixture와 전체 결과는
`docs/report/work/complex-case`와 `docs/report/work/complex-run`에 있다. 보고서에 사용하는
근거는 `docs/report/evidence/complex-*`와 `benchmark-*.json`에 보존했다.

| 입력 | 전극 / 전극 조합 | TAT 중앙값 | 최대 RSS | 출력 크기 |
|---|---:|---:|---:|---:|
| pilot 200×200 | 4 / 3 | 1.24초 | 196.4 MB | 0.08 MB |
| complex 200×200 | 10 / 16 | 2.24초 | 218.5 MB | 14.00 MB |
| pilot 1000×1000 | 4 / 3 | 3.36초 | 420.5 MB | 0.95 MB |
| complex 1000×1000 | 10 / 16 | 13.48초 | 745.1 MB | 65.16 MB |

- 1000×1000에서 복합/pilot TAT는 4.01배, 최대 RSS는 1.77배였다.
- 복합 입력을 200×200에서 1000×1000으로 키우면 TAT는 6.02배였다.
- 각 수치는 1회 워밍업 뒤 5회 실행한 결과이며 TAT는 중앙값, RSS는 최대값이다.
- 지정한 vertical short 양성 전극 조합 8개는 모두 0이 아니었고, 먼 거리 음성 대조 2개는 모두 0 px였다.
- 모든 이물 유입 공정·전극 조합에서 크기 단조성 위반 0, 중복 제거/동시 영향 조합 수 불일치 0, CA 적분 오차 0이었다.
- 한 픽셀에서 최대 4개 전극 조합의 hotspot이 겹쳐 중복 제거 경로가 실제로 사용됐다.
- 내부 기존 방식/EDT 차이 2,911 px는 kernel 모양 규칙으로 모두 분류됐고 미분류는 0이었다.

기존 현장 `sim` 실행 파일과 동일 입력 변환기는 이 저장소에서 찾지 못했으므로 외부 기존
산출물과의 픽셀 대조는 완료되지 않았다. 현장 실행 환경을 확보하면 단계 C를 이어서 수행한다.

2026-09-19에는 같은 200×200 복합 형상에 open 경로 3개를 추가한 `2B` 파일럿도 수행했다.
5회 측정 중앙값은 lateral short 1.444초, vertical short 1.785초, open 16.523초, 통합 17.862초다.
단독 실행 TAT 합에서 open이 83.65%를 차지했으며, 통합 CA 구성비는 lateral short 3.48%,
vertical short 61.26%, open 35.27%였다. 이 CA 비율은 Fab 가중 전 합성 기하 결과다.

같은 형상·크기의 open 결과를 적용 공정 사이에서 재사용한 뒤 open은 7.545초, 통합은 8.810초로
각각 54.3%, 50.7% 감소했다. 최적화 전후 hotspot 전체 배열과 CA 값은 동일했다. 수정 후에도
open이 단독 TAT 합의 70.25%였으므로, 두 번째 고속화는 후보마다 전체 ROI 연결 성분을 다시 계산하는
부분을 대상으로 했다.

연결 성분 판정을 전극 bounding box로 제한한 최종 결과는 open 3.135초, 통합 4.423초다. 최초
기준선 대비 각각 5.27배와 4.04배 빨라졌으며, open의 단독 TAT 비율은 49.43%로 낮아졌다. 최초
기준선과 최종 결과의 hotspot·CA는 완전히 같았다.

같은 국부 형상을 1000×1000 중앙에 둔 확장성 시험에서는 lateral short 5.678초, vertical short
8.522초, open 4.669초, 통합 16.385초였다. 통합 최대 RSS는 843.6MB, 출력은 621.4MB였다.
200×200 대비 통합 TAT는 3.70배였고 open은 1.49배였다. 절대 CA는 동일했으며, 출력의 대부분은
float32 `pair_margin_um` 지도였다. 이 결과는 ROI와 출력 크기 증가를 보여주지만 실제 긴 배선의
open TAT를 대신하지 않는다.
