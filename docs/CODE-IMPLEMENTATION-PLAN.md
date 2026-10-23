# particle-sim 코드 보완 실행 계획

작성일: 2026-09-19  
기준: 지금까지 합의한 해석 범위, 현재 소스 코드, 합성 2nd test와 보고서 기록

이 문서는 `docs/CODE-CHANGE-BACKLOG.md`의 상세 근거를 실제 구현 순서로 다시 배열한 계획이다.

## 구현 상태 · 2026-09-19

P0 코드 보완은 완료했다. layer 잔여 영역, 전극 조합 전수 감사, 단자 영역 open과 적용 공정,
x/y overlay·항목별 근거, 명목/산포 CA 분리, `1/k` 기여도, Fab 평균 이물 입력과 네 가지 위험도,
artifact hash, 대형 산포 mask 저장 방식, 외부 기존 sim 비교기를 구현하고 자동 검증을 통과했다.

1000×1000, 64조건 기존 성능 fixture 재측정은 wall 107.325초, 최대 RSS 636.8MB였다. 이전
102.254초/1.072GB 대비 TAT는 약 5% 늘고 최대 RSS는 약 41% 줄었다. 기능 계약의 P0 구현은 끝났지만,
실제 기존 sim 산출물이 저장소에 없어 **실제 동일 입력 비교와 대체 판정은 아직 남아 있다**. 자료가
확보되면 이 비교를 최우선으로 수행하고, 그 결과가 확정된 뒤 simcompare 계약을 연결한다. 보고서
HTML/PDF는 마지막에 한 번 재생성한다.

복합 200×200 형상의 failure mode별 파일럿도 완료했다. 1회 워밍업 후 5회 중앙값은 lateral short
1.444초, vertical short 1.785초, open 16.523초, 통합 17.862초였다. open은 활성 계산 조건이
35개로 가장 적어도 단독 TAT 합의 83.65%를 차지했다. 동일 형상·크기의 open 결과를 적용 공정
사이에서 재사용한 뒤 open은 7.545초, 통합은 8.810초로 각각 54.3%, 50.7% 감소했다. 최적화 전후
hotspot·CA 픽셀 동등성과 Zarr hash 일치를 확인했다. 수정 후에도 open 비중이 70.25%이므로
연결 성분 판정을 전극 bounding box로 제한했다. 최종 open은 3.135초, 통합은 4.423초로 최초
기준선보다 각각 5.27배와 4.04배 빨라졌고 픽셀 동등성을 다시 확인했다.

같은 200×200 형상을 1000×1000 중앙에 그대로 둔 확장성 시험도 5회 수행했다. lateral short
5.678초, vertical short 8.522초, open 4.669초, 통합 16.385초였고 통합 최대 RSS는 843.6MB였다.
200×200 대비 통합은 3.70배, open은 1.49배였다. open 연결성은 전극 bounding box만 계산하므로
ROI 빈 공간 증가의 영향이 제한됐지만, short 지도는 전체 ROI에 저장되므로 통합 출력이 621.4MB로
늘었다. 그중 `pair_margin_um`이 621.3MB(592.5MiB)로 대부분을 차지했다. 이는 같은 국부 형상을 큰
canvas에 둔 ROI·출력 확장성 시험이며, 긴 배선이나 실제 제품 형상을 대표하지 않는다.

### 2026-09-19 검토 후 추가 계획

미커밋 작업 트리 검토에서 open 판정 결함 두 가지를 확인했다. 입자가 없어도 이미 끊긴 전극, 그리고
산포 변형으로 단자가 사라진 조건에서 후보 위치 전부가 open으로 계산된다. 이 수정과 함께 open 단자를
contact·via mask에서 자동으로 만들고, contact가 이물에 막히는 **CNO(Contact Not Open)**를 open,
short와 별도 불량 유형으로 추가한다(2.3절, 단계 5a). 실제 기존 sim 자료 확보와 병행하며, 실제 제품
mask 입력 전에 완료한다. 나머지 검토 결함은 `CODE-CHANGE-BACKLOG.md` 14절에 있다.

**단계 5a 구현 상태(2026-09-19):** 1~7번을 구현했다. 전체 자동 시험은 297개가 통과했다(신규 18개 포함).
contact가 없는 기존 복합 200² 파일럿은 hotspot Zarr 5개 배열과 CA·union 표가 이전 결과와 모두 같고,
`ca.parquet`에 `contact_id` 열만 추가됐다. 같은 파일럿의 G1·D1·D2 단자 영역을 같은 net의 pad 전극과
contact hole로 바꾼 사례에서 자동 단자 open은 명시 단자 결과와 같았다. 차이는 contact가 완전히 막힌
위치(G1 552 px, D1 368 px, D2 368 px)뿐이며 모두 CNO로 옮겨 갔다. 5회 중앙값은 명시 단자 3.724초,
contact 단자 3.889초(+4.4%), 최대 RSS는 224.8MB와 243.5MB다. 이 측정은 프로세스 안의 `run_nominal`
시간이라 CLI 전체를 잰 기존 표(4.423초)와 직접 비교하지 않는다. 근거:
`docs/report/evidence/benchmark-contact-pilot-200.json`, 재실행 `docs/report/benchmark_contact_pilot.py`.

같은 두 사례에 가정 산포(비기준 layer overlay x/y σ 0.2 µm, CD σ 0.1 µm)로 Sobol 64조건을 1회씩
실행했다. 명시 단자 325.3초·최대 RSS 330.3MB, contact 단자 287.8초·377.5MB였다. 1회 측정이라
두 값의 차이를 성능 차이로 해석하지 않는다. 두 사례 모두 산포 제외 조건(`baseline_open`,
`terminal_lost`, `contact_lost`)과 소실 contact hole은 0건이었다. contact 사례의 CNO 평균 CA는
177.3 µm²(명목 138.0 µm²), 자동 open 평균 CA는 6,088.4 µm²(명목 6,071.0 µm²)다. 근거:
`docs/report/evidence/benchmark-contact-probability-200.json`, 재실행
`docs/report/benchmark_contact_probability.py`.

2026-09-20에 검토 14절의 나머지 결함(비교 계약 식별자, 부분 출력, 수렴 지표, 평균 이물 표 경로,
임시 폴더)도 수정했고 전체 자동 시험 304개가 통과했다. 명목에서 이미 단자가 부족한 단위는 산포
제외 건수에 세지 않는다. 실제 데이터 시험 전 단계의 코드 작업은 모두 끝났다.

## 1. 현재 상태 판정

기하 계산의 기본 방향은 유지할 수 있다.

- 현재 sim도 높이를 쓰는 2.5D 계산이다.
- 개선 핵심은 기존 계산 대상에 없던 전극 조합과 xy 거리와 z 거리를 함께 쓰는 대각 CA다.
- `lateral short`, `vertical short`, `open`, `hotspot` 용어는 그대로 유지한다.
- hotspot은 이물 중심 후보 위치이며, 산포 계산에서는 각 위치의 조건부 발생 비율로 해석한다.
- `union:any`는 한 layer의 전극 형상 합집합이 아니다. 같은 이물 유입 공정과 크기에서 선언된
  모든 `lateral short`, `vertical short`, `open` hotspot을 픽셀별 OR한 결과다.

현 상태로 실제 제품과 Fab 데이터를 바로 연결하면 안 된다. 아래 문제가 결과 의미를 바꿀 수 있다.

| 구분 | 보완 전 동작 | 영향 | 상태 |
|---|---|---|---|
| layer와 전극 mask | 전극 합집합과 layer mask의 완전 일치를 강제 | TR/기타 패턴이 남는 실제 mask를 거부 | 완료 |
| 전극 조합 | YAML에 쓴 조합만 계산 | 누락 조합을 hotspot 0과 구분할 수 없음 | 완료 |
| open 단자·공정 | 점 단자와 모든 공정 반복 | 실제 입력과 공정별 위험이 모호 | 완료 |
| overlay·근거 | layer당 동일 x/y와 파일 단위 근거 | 실측 산포를 구분하지 못함 | 완료 |
| union CA·Fab | 명목/평균 혼용과 pair 중복 집계 | 위험도 합계와 의미가 틀릴 수 있음 | 완료 |
| 위험도 출력 | layer·전극 조합·전극 기여도가 없음 | 주요 전극을 찾는 목적을 충족하지 못함 | 완료 |
| 실행 이력·메모리 | 출력 hash 없음, 모든 변형 mask 선보관 | 재현성과 대형 입력 메모리 문제 | 완료 |
| 기존 sim 비교 | 외부 결과 importer와 픽셀 비교 명령이 없음 | 대체 여부를 판단할 수 없음 | 엔진 완료, 실제 자료 대기 |
| simcompare | 실제 전극 ID와 조합 ID를 구조화해 전달하지 않음 | 제품 내부 전극 순위가 유지되지 않음 | P1 대기 |

보완 전 회귀 기준은 `pytest -q` **257 passed**였고, 보완 후 기준은 **279 passed**다. 기존 측정
기록은 다음과 같다.

| 실행 | 중앙 TAT | 최대 RSS | 비고 |
|---|---:|---:|---|
| pilot 200² 기준 조건 | 1.238초 | 196MB | 5회 측정 |
| 복합 200² 기준 조건 | 2.237초 | 219MB | 전극 10개, 조합 16개 |
| pilot 1000² 기준 조건 | 3.361초 | 421MB | 5회 측정 |
| 1000² 산포 64조건 | 102.254초 | 1.07GB | 보완 전 eager mask 생성 포함 |
| 복합 200² lateral short | 1.444초 | 198MB | 6조합, 활성 조건 90개, 5회 측정 |
| 복합 200² vertical short | 1.785초 | 206MB | 10조합, 활성 조건 150개, 5회 측정 |
| 복합 200² open · 기준선 | 16.523초 | 194MB | 3경로, 활성 조건 35개, 5회 측정 |
| 복합 200² 통합 · 기준선 | 17.862초 | 225MB | short 16조합+open 3경로, 5회 측정 |
| 복합 200² open · 공정 간 cache | 7.545초 | 194MB | 기준선 대비 54.3% 감소, 5회 측정 |
| 복합 200² 통합 · 공정 간 cache | 8.810초 | 222MB | 기준선 대비 50.7% 감소, 5회 측정 |
| 복합 200² open · 전극 bbox | 3.135초 | 194MB | 최초 기준선 대비 81.0% 감소, 5회 측정 |
| 복합 200² 통합 · 전극 bbox | 4.423초 | 222MB | 최초 기준선 대비 75.2% 감소, 5회 측정 |
| 동일 형상 padded 1000² lateral short | 5.678초 | 534MB | ROI·출력 확장성 시험, 5회 측정 |
| 동일 형상 padded 1000² vertical short | 8.522초 | 650MB | ROI·출력 확장성 시험, 5회 측정 |
| 동일 형상 padded 1000² open | 4.669초 | 333MB | 전극 bbox 유지, 5회 측정 |
| 동일 형상 padded 1000² 통합 | 16.385초 | 844MB | 출력 621.4MB, 5회 측정 |

복합 200²는 단순 pilot보다 약 81% 오래 걸렸다. 이후 변경은 같은 입력과 측정 방법으로 TAT,
최대 RSS, 출력 크기를 다시 측정한다.

## 2. 이번 구현에서 고정할 계산 의미

### 2.1 입력 mask

- 각 전극 mask는 대응 layer mask의 부분집합이어야 한다.
- 같은 layer의 서로 다른 전극 mask는 겹치면 안 된다.
- `layer mask - 전극 mask 합집합`은 허용하며 픽셀 수, 면적, layer 대비 비율을 기록한다.
- 남는 영역은 현재 계산에서 `TR 또는 기타 미평가 영역`이다. TR의 신호 분리 기능은 아직
  `lateral short`, `vertical short`, `open` 판정에 넣지 않는다.
- 표시 예제에서는 어려운 약어 대신 `mask1`, `mask2`를 사용한다. 실제 layer ID는 입력값을 보존한다.

### 2.2 전극 조합 누락 검사

모든 다른 net 전극 조합을 먼저 목록화한다. 같은 layer면 `lateral short`, 다른 layer면
`vertical short` 후보다. 기하학적 거리가 멀다는 이유만으로 자동 제외하지 않는다. 이물 크기,
overlay와 CD 산포가 바뀌면 도달 가능성이 달라질 수 있기 때문이다.

각 후보는 다음 중 하나여야 한다.

1. `electrode_pairs`에 계산 대상으로 선언
2. 제외 규칙 또는 제외 목록에 사유와 함께 선언

`validate-input`은 선언, 제외, 미결정 조합을 표와 JSON으로 내고 미결정 조합이 있으면 실제 실행을
막는다. 기하학적 최소 xy 거리와 layer 높이 차이는 검토용 수치로만 함께 출력한다. 먼 거리 음성
대조는 계산 대상으로 계속 선언할 수 있다.

### 2.3 open

실제 제품은 점 좌표 대신 두 단자 영역 mask를 기본 입력으로 쓴다.

- 단자 영역은 비어 있지 않고 대상 전극 영역과 교차해야 하며 서로 겹치면 안 된다.
- layer overlay는 단자 영역에도 같은 이동을 적용한다. CD는 전극 형상에 적용하고 단자 영역 자체를
  팽창·수축하지 않는다.
- 이물이 가린 뒤 각 단자 영역에 남은 전극 픽셀 중 하나라도 서로 연결돼 있으면 정상이다.
- 어느 한 단자 영역에 남은 전극 픽셀이 없거나 두 영역을 잇는 경로가 없으면 `open`이다.
- 기존 단일 픽셀 단자는 합성 fixture와 v1 입력 호환용으로만 유지한다.
- `roi_boundary`는 `open proxy`로 표시하며 실제 제품 확정 결과에는 사용할 수 없게 한다.

또한 각 open 경로에 적용 가능한 이물 유입 공정을 명시한다. 현재 코드는 source와 무관한 같은 open을
모든 공정에 반복하므로, 실제 입력에서는 `applicable_particle_source_steps`가 없으면 실패시킨다.
합성 v1 입력만 호환 경고와 함께 기존 동작을 허용한다.

#### 2.3.1 contact 기반 자동 단자 open (추가 계획)

- `contacts` 입력으로 contact·via 원본 mask와 `hole_layer`를 받는다. 위·아래 도체 layer는 stack
  `process_order`로 정하고, 여러 절연층을 관통하면 직접 지정한다. 이어지는 두 전극은 같은 net이어야 한다.
- 전극별 자동 단자는 겹치는 contact hole 하나하나와 ROI 경계에 닿은 구간이다. ROI 경계는 배선이 ROI
  밖으로 이어지므로 단자로 본다. 단자가 2개 미만이면 `insufficient_terminals`로 기록한다.
- 입자가 없는 상태에서 단자를 연결 성분별로 묶은 기준 분할을 먼저 만든다. 이물로 같은 묶음이던
  단자가 분리되거나 ROI 경계 단자가 완전히 덮이면 open이다. 단자 없는 조각의 분리는 무시한다.
- contact 단자가 막히는 경우는 open이 아니라 CNO로 계산하고, 그 위치의 open 비교에서 뺀다.
- 판정은 전극 단위다. 같은 net의 다른 layer 우회 경로를 반영하는 net 단위 판정은 후속 항목이다.
- 적용 이물 유입 공정은 layer별 `open_applicable_steps_by_layer`로 명시한다. 적층 순서에서 자동으로
  추정하지 않는다.
- 명시한 `open_paths`는 자동 단자보다 우선한다. 명시 경로에도 입자 없는 상태의 연결 검사를 적용한다.
  명목 입력에서 끊겨 있으면 입력 오류, 산포 조건에서 끊기면 `baseline_open`으로 기록한다.
- 산포에서 contact는 `hole_layer`의 overlay·CD를 따른다. 변형 후 단자와 전극의 교집합이 비면
  `terminal_lost`로 기록하고 open으로 세지 않는다.

#### 2.3.2 CNO (추가 계획)

- CNO는 이물이 contact·via hole을 막아 위·아래 전극이 연결되지 않는 불량이다.
- 불량 단위는 contact 연결 하나(`cno:<contact_id>:<lower>__<upper>`), `failure_mode=contact_not_open`,
  `unit_kind=contact`다.
- 잠정 기준: hole이 **전부** 가려질 때 CNO, 같은 연결의 hole 중 하나라도 막히면 CNO. 일부 가림과
  병렬 hole 기준은 논의가 필요한 개선 항목이며, 실행 manifest와 보고서에 잠정 기준임을 표시한다.
- 적용 이물 유입 공정은 contact별로 명시한다. 추천 선언 범위는 아래 도체 형성 후부터 위 도체
  증착 전까지다.
- 산포에서 hole이 전극과 겹치지 않게 되면 `contact_lost`로 기록하고 CNO CA에 넣지 않는다.
- union에 `contact_not_open` 유형을 추가하고 `union:any`에 포함한다. 전극 기여도는 위·아래 전극에
  절반씩 배분하고 `risk_by_contact.parquet`를 추가한다.

### 2.4 산포와 Sobol

- 첫 실제 산포 입력은 Fab 실측 mask layer별 overlay x, overlay y, CD만 사용한다.
- overlay x/y, mask별 overlay/CD는 현재 확인대로 서로 독립으로 둔다.
- 기준 layer, 상대 이동 부호, CD가 전체 선폭 변화인지 한쪽 경계 변화인지 입력에 기록한다.
- 각 분포마다 `measured` 또는 `assumed`, 원자료 ID/hash, 단위, 제품군, 공정, lot/기간을 기록한다.
- 두께, 이물 직경, 이물 중심 offset, 매립 깊이, 평탄화와 접촉 여유는 실측 근거가 생길 때까지
  0 또는 기준값으로 고정하고 실측 산포처럼 표시하지 않는다.
- 첫 비교의 `contact_margin_um`은 0으로 고정한다.

Sobol은 정규분포를 대신하지 않는다. Sobol 점으로 0~1 공간을 고르게 채운 뒤 각 점을 정규분포 또는
실측 경험분포의 역누적분포에 통과시킨다. 따라서 최종 표본의 분포는 입력 분포를 따른다.

수렴은 64→128→256 조건을 기본으로 보고, 서로 다른 scramble seed를 최소 3개 사용한다. 256에서
기준을 만족하지 않으면 512로 올린다. 주 판정값은 각 이물 유입 공정·크기의 `union:any`이며,
전극 기여도 상위 항목의 순위 안정성도 함께 확인한다. 허용 오차는 pilot 결과와 함께 명시적으로
고정하고 코드에 숨은 기본값을 두지 않는다.

### 2.5 중복 없는 전극 기여도

산포 조건 `s`, 실패 단위 `i`, 픽셀 `x`에서 hotspot을 `H(i,s,x)`라 한다. 같은 픽셀의 활성 단위
개수 `k(s,x)`를 구하고 다음 값을 산포 조건마다 먼저 계산한다.

- 원본 민감도: `H(i,s,x)`의 전체 면적. 다른 조합과 합산할 수 없다.
- 배분 기여도: hotspot 픽셀마다 `1 / k(s,x)`를 해당 활성 단위에 배분한다.
- 고유 기여도: `k(s,x) == 1`이면서 해당 단위만 활성인 면적이다.

모든 실패 단위의 배분 기여도 합은 같은 조건의 `union:any` CA와 정확히 같아야 한다. 산포 평균과
분위수는 이 배분을 각 조건에서 수행한 뒤 계산한다.

`open`을 전극 조합이라고 부르지 않도록 내부 기준 ID는 `failure_unit_id`로 통일한다.

- `risk_by_pair.parquet`: `lateral short`, `vertical short`만 포함
- `risk_by_open_path.parquet`: `open`만 포함
- `risk_by_electrode.parquet`: short는 조합 배분 기여도의 절반씩 두 전극에, open은 대상 전극에 전부 배분
- `risk_by_layer.parquet`: 이물 유입/Fab layer 기준 전체 union과 영향을 받는 전극 layer 기준 배분값

전극 및 영향 layer 배분은 물리적 인과가 아니라 합계 중복을 없애기 위한 보고 규칙임을 metadata에
기록한다. 원본 민감도, 배분 기여도, 고유 기여도를 모두 남겨 순위 해석을 한 값에 의존하지 않게 한다.

### 2.6 Fab 이물 자료

보유 데이터는 layer별·이물 크기 구간별 평균 이물 개수이며 검사면적도 알고 있다. 첫 지원 계약은
다음 값으로 고정한다.

- `fab_distribution_layer`
- `size_bin_left_um`, `size_bin_right_um`, `size_representative_um`
- `mean_particle_count`
- `inspection_area_per_sample_mm2`
- `mean_density_per_mm2 = mean_particle_count / inspection_area_per_sample_mm2`
- `sample_count` 선택값, `measurement_period`, `particle_type`, 원자료 ID/hash

`source_step` 문자열로 Fab 표를 직접 찾지 않고 manifest의 명시적인 `fab_distribution_layer`로 찾는다.
한 크기 bin의 평균 개수를 여러 시뮬레이션 직경에 그대로 반복 적용하지 않는다. 첫 구현은 한 bin과
한 대표 직경이 정확히 대응하는 방식만 실제 Fab 모드에서 허용한다. 보간이나 bin 내 적분은 밀도 질량을
보존하는 규칙을 별도 검증한 뒤 추가한다.

가중값은 `CA(mm²) × mean_density(1/mm²)`인 기하학적 hit 기대값이다. 실제 불량률이나 수율로 부르지
않는다. Fab 적용은 크기별 기하 시뮬레이션이 끝난 뒤 수행하며 Sobol의 이물 크기 산포와 섞지 않는다.

## 3. 출력 schema 2.0

입력 v1은 명시적인 변환 경고와 함께 내부 v2 모델로 읽을 수 있게 한다. 새 실행은 schema 2.0 출력만
생성한다. 기존 출력 파일을 조용히 새 의미로 읽지 않고 변환 명령 또는 이해하기 쉬운 오류를 낸다.

### 기본 기하/산포 출력

- `ca.parquet`
  - `failure_unit_id`, `unit_kind`, `electrode_pair_id`, `open_path_id`
  - 실제 전극 ID, 전극 layer, net, failure mode
  - 기준 CA와 산포 평균·p05·p50·p95를 별도 열로 보존
- `ca_union.parquet`
  - `union_scope=all_failure_modes|failure_mode`
  - `ca_union_nominal_um2`
  - `ca_union_mean_um2`, `ca_union_p05_um2`, `ca_union_p50_um2`, `ca_union_p95_um2`
  - `scenario_count`, `convergence_delta`
  - 단순 합은 진단용 `ca_raw_unit_sum_*`으로만 보존
- `hotspots.zarr` / `probability.zarr`
  - 축, 단위, `failure_unit_id`, mode, 전극 ID와 union 범위를 attrs에 기록
  - Zarr는 큰 다차원 배열을 chunk 단위로 읽고 쓰는 계산용 저장 형식임을 사용자 문서에 설명
- Parquet 표
  - 열 이름·단위·자료형을 보존하고 필요한 열만 빠르게 읽는 분석용 표 형식임을 사용자 문서에 설명

산포가 0이면 nominal과 mean이 같아야 한다. 산포 실행에서 mean을 nominal 열에 쓰는 현재 동작은
schema 2.0에서 제거한다.

### Fab 및 순위 출력

| 파일 | 주 용도 | 합산 가능 값 |
|---|---|---|
| `risk_by_layer.parquet` | 이물 유입 layer 및 영향 layer별 위험 | union/배분 기여도 |
| `risk_by_pair.parquet` | short 전극 조합 순위 | 배분 기여도 |
| `risk_by_open_path.parquet` | open 경로 순위 | 배분 기여도 |
| `risk_by_electrode.parquet` | 개별 전극 순위 | allocated electrode contribution |
| `risk_summary.csv` | 사람이 바로 확인하는 상위 순위와 단위 | 위 결과의 요약 |

각 상세 행은 이물 유입 공정, Fab layer, 이물 크기, 실제 전극 ID와 layer, failure mode를 함께 가진다.
따라서 layer 위험, 전극 위험, layer×전극 위험을 같은 실행에서 추적할 수 있다.

### 실행 이력

`run_manifest.json`에 다음을 추가한다.

- git commit과 working tree dirty 여부
- git을 쓸 수 없을 때 package/build ID
- 입력 schema, 출력 schema와 변환 이력
- 실행 명령과 핵심 옵션, Sobol seed 목록
- 입력 파일 hash
- 주요 Zarr/Parquet/CSV의 상대 경로, 크기와 hash
- 산포 항목별 evidence와 원자료 식별자
- TR 미평가, contact margin 0, 원형/구형 이물 등 현재 제한

manifest 자신의 hash는 artifacts 목록에서 제외해 순환 참조를 피한다. Zarr는 정렬한 상대 파일명과
각 파일 내용으로 디렉터리 digest를 만든다.

## 4. 구현 순서

### 단계 0 — 현재 기준 보존

대상: 현재 fixture, 테스트, benchmark 기록

1. 보완 전 257개 테스트 결과와 합성 출력 schema를 기준 artifact로 보존하고, 보완 후 279개 회귀 결과를 추가한다.
2. pilot 200², 복합 200², pilot 1000², 1000² 산포 64조건 측정 명령을 재사용 가능한 benchmark
   script로 옮긴다.
3. 기존 `legacy_failure`는 외부 기존 sim 결과가 아니라 내부 형태연산 기준임을 이름과 문서에서
   구분한다.

완료 조건: 동일 commit과 입력 hash로 기준 수치와 파일을 다시 만들 수 있다.

### 단계 1 — schema와 입력 검증

대상 파일:

- `src/particle_sim/schema.py`
- `src/particle_sim/config.py`
- `src/particle_sim/io/masks.py`
- 신규 `src/particle_sim/validation/electrode_pairs.py`
- `src/particle_sim/cli.py`
- `SPEC.md`, `ARCHITECTURE.md`, `docs/arch/*`

작업:

1. layer–전극 관계를 완전 일치에서 부분집합 검사로 바꾸고 잔여 영역 통계를 저장한다.
2. 모든 다른 net 전극 조합을 목록화하고 선언/제외/미결정 상태를 검사한다.
3. open 단자 영역 mask와 적용 가능한 이물 유입 공정 목록을 추가한다.
4. overlay x/y 분리, 기준 layer, 독립 관계와 항목별 evidence/provenance를 추가한다.
5. 사용자 표시에서는 `source_step`, `pair` 같은 축약 대신 `이물 유입 공정`, `전극 조합`을 쓴다.

완료 조건:

- TR 잔여 영역은 통과하고 layer 밖 전극은 실패한다.
- 미결정 전극 조합이 하나라도 있으면 실제 실행 전에 실패한다.
- 단자 mask의 shape, 공집합, 전극 교차, 상호 중첩 오류를 검출한다.
- v1 fixture 결과는 호환 경로에서 픽셀 단위로 유지된다.

### 단계 2 — 산포·union·기여도 계산 정정

대상 파일:

- `src/particle_sim/uncertainty/sampler.py`
- `src/particle_sim/geometry/transforms.py`
- `src/particle_sim/aggregation/probability.py`
- 신규 `src/particle_sim/aggregation/risk_allocation.py`
- `src/particle_sim/failure/union.py`
- `src/particle_sim/units.py`
- `src/particle_sim/pipeline.py`
- `src/particle_sim/io/outputs.py`
- `src/particle_sim/validation/convergence.py`

작업:

1. mask별 overlay x/y와 CD를 각각 독립 Sobol 차원에서 역누적분포 변환한다.
2. nominal union과 산포 mean/분위수 열을 분리한다.
3. 각 조건의 hotspot OR, `1/k` 배분, 고유 기여도를 한 번의 픽셀 순회로 계산한다.
4. short pair와 open path를 공통 `failure_unit_id`로 처리하되 출력에서는 분리한다.
5. open을 지정된 이물 유입 공정에서만 계산하고 안전한 조건에만 캐시한다.
6. 3개 이상의 scramble seed와 64→128→256→선택 512 수렴 기록을 지원한다.
7. 모든 조건의 변형 mask를 미리 메모리에 보관하지 않도록 순차 처리 또는 제한된 cache로 바꾼다.

완료 조건:

- `sum(allocated failure-unit CA) == union:any CA`가 각 조건과 평균에서 성립한다.
- 산포 0에서 nominal과 mean이 일치한다.
- x/y 분포 표본 통계와 독립성이 입력값에 맞는다.
- 1000² 산포 64조건의 최대 RSS가 기존 1.07GB보다 감소하고 hotspot/CA는 변경 이유가 없는 항목에서
  동일하다.

### 단계 3 — Fab 연결과 위험도 출력

대상 파일:

- `src/particle_sim/io/fab.py`
- `src/particle_sim/aggregation/fab_weighting.py`
- 신규 `src/particle_sim/aggregation/risk_tables.py`
- `src/particle_sim/cli.py`
- `src/particle_sim/visualization/reports.py`

작업:

1. Fab 입력을 평균 이물 개수/샘플 검사면적 계약으로 바꾼다.
2. `particle_source_step → fab_distribution_layer` 매핑을 사용한다.
3. 한 bin을 한 대표 직경에만 연결하고 밀도 중복 사용을 막는다.
4. union 전체값, pair/open 원본·배분·고유 기여도, 개별 전극과 영향 layer 배분값을 만든다.
5. `apply-fab --run <실행 폴더>`가 필요한 CA·union·manifest를 함께 읽고 다섯 위험도 파일을 한
   출력 폴더에 생성하도록 바꾼다.

완료 조건:

- 평균 이물 개수/검사면적 수기 계산과 밀도가 일치한다.
- `source_step` 이름과 Fab layer 이름이 다른 fixture에서도 지정한 Fab 자료만 사용한다.
- Fab 가중 후 layer별 union 합과 pair/open 배분 합, 전극 배분 합이 각각 일치한다.
- 원본 민감도는 합산 불가라고 schema metadata와 요약표에 표시한다.

### 단계 4 — 재현성과 schema 고정

대상 파일:

- `src/particle_sim/io/outputs.py`
- 신규 provenance/hash helper
- 모든 출력 reader와 report

작업:

1. schema 2.0 writer/reader와 v1 입력 변환 경로를 완성한다.
2. commit, dirty 상태, build ID와 artifact hash를 기록한다.
3. Zarr attrs, Parquet 열과 manifest 간 ID·단위 일치 검사를 추가한다.
4. 출력 schema 문서와 예제 YAML/CSV를 실제 생성 결과에서 자동 검증한다.

완료 조건: 결과 폴더만 전달해도 코드 버전, 입력, 출력 무결성, 계산 의미를 확인할 수 있다.

### 단계 5 — 기존 sim 비교기 구현과 로컬 비교

대상 파일:

- 신규 `src/particle_sim/comparison/legacy_sim.py`
- `src/particle_sim/cli.py`
- 비교 fixture와 결과 report

먼저 기존 sim 결과를 다음 공통 형식으로 가져온다.

- 배열 축과 원점, pixel 크기
- 이물 크기가 직경인지 반경인지
- 이물 유입 공정, 실제 전극 조합 ID, failure mode
- hotspot 배열 또는 CA, 경계 처리와 kernel 정보

비교 출력은 전극 조합·이물 크기별 XOR map, 교집합, 기존 sim 전용/현재 sim 전용 픽셀, CA 절대·상대
차이와 원인 분류를 포함한다. 두 sim이 공통으로 계산한 조합만 증감 비교한다. 기존 계산 대상에 없던
대각 전극 조합은 기존값 0이 아니라 `N/A`, `new_pair`로 표시한다. 차이는 “기존 2D에서 2.5D로
변경”으로 설명하지 않고 “기존 2.5D에 없던 xy+z 대각거리와 신규 조합 추가”로 분류한다.

완료 조건:

- 공통 형식 round-trip에서 XOR가 0이다.
- 축, 단위, 크기 정의, 누락 조합을 자동 추정하지 않고 오류로 낸다.
- 실제 동일 입력 비교에서 미분류 차이가 남으면 대체 판정을 통과시키지 않는다.
- 비교한 실행은 clean commit 또는 명확한 dirty 상태와 artifact hash를 가진다.

기존 sim의 실제 파일 형식 adapter를 완성하려면 대표 입력과 출력 한 세트, 축/단위/직경 정의가
필요하다. 이 자료가 준비되기 전에도 공통 importer와 비교 엔진까지는 구현할 수 있다.

### 단계 5a — open 판정 정정, 자동 단자 open, CNO

대상 파일: `schema.py`, `config.py`, `io/masks.py`, `geometry/transforms.py`,
`failure/open_failure.py`, 신규 `failure/contact_not_open.py`, `units.py`, `failure/union.py`,
`aggregation/risk_allocation.py`, `aggregation/risk_tables.py`, `io/outputs.py`, `pipeline.py`,
`comparison/legacy_sim.py`, 관련 테스트와 입력·출력 가이드

1. 명시 `open_paths`에 입자 없는 상태의 연결 검사와 `baseline_open`·`terminal_lost` 분류를 추가한다.
   검토 결함 1·2의 재현 사례를 회귀 테스트로 고정한다.
2. `contacts` schema와 입력 검증(위·아래 전극 겹침, 같은 net, hole 단위 분리)을 추가한다.
3. 자동 단자(contact hole, ROI 경계)와 기준 분할 기반 전극 단위 open을 구현한다.
4. CNO 판정, 적용 공정, 산포의 `contact_lost` 처리를 구현한다.
5. union·`1/k` 배분·Fab 위험도·출력 schema에 CNO와 자동 open 단위를 추가한다.
6. 기존 sim 비교에서 새 단위를 `new_unit`으로 분리한다.
7. 복합 200² 파일럿에 contact mask를 추가해 5회 TAT·RSS를 측정하고 기존 기록과 비교한다.

완료 조건: 명시 단자 영역을 contact로 바꾼 자동 단자 결과가 기존 방식과 픽셀 일치한다(막힌 contact
위치가 CNO로 옮겨 간 차이만 허용). 기준 단선·단자 소실 사례에서 후보 전체가 open이 되지 않는다.
open·short·CNO가 겹치는 사례에서 배분 합계와 `union:any`가 일치한다.

### 단계 6 — 실제 제품 pilot와 simcompare 연결

기존 sim 비교가 끝난 뒤 다음 순서로 진행한다.

1. 실제 제품 nominal 계산: `contact_margin_um=0`, 원형/구형 이물, TR 전기 기능 미반영
2. 실제 mask별 overlay/CD 산포 적용: 3 seeds, 64→128→256, 필요 시 512
3. 크기별 시뮬레이션 완료 후 Fab 평균 이물 개수 적용
4. layer, 전극 조합, open 경로, 개별 전극, layer×전극 순위 확인
5. `particle-sim-compare` 계약에 실제 전극 ID, 전극 조합 ID, 이물 유입 공정, Fab layer와
   산포/기준 hotspot 의미를 전달

simcompare의 실측 평균 이물 개수 모드와 전극 ID별 순위가 구현되기 전까지 exporter 성공만으로
연동 완료라고 판단하지 않는다.

## 5. 검증 묶음

필수 자동 검증은 구현과 함께 추가한다.

1. mask 부분집합, 잔여 TR 통계, layer 밖 전극 오류
2. 전극 조합 선언/제외/미결정 전수 검사와 순서 독립성
3. 단자 영역 일부 가림, 전체 가림, 우회 경로, 분기, ㄱ자 open 전수 비교
4. source별 open 적용 여부
5. overlay x/y 역누적분포, 독립 차원, 기준 layer와 부호
6. nominal/mean 열 분리와 산포 0 일치
7. 겹침 없음, 완전 겹침, 부분 겹침에서 `1/k` 배분 불변식
8. open 포함 전극 및 layer 배분 합계
9. 평균 이물 개수/면적 밀도와 대표 직경 1:1 연결
10. 출력 schema, ID, 단위, artifact hash와 v1 호환 오류
11. 기존 sim 공통 형식 round-trip과 축/단위 오류
12. 기존 257개 회귀 테스트와 신규 계약 테스트
13. 기준 단선(두 섬 전극), 산포 단선, 단자 소실에서 `baseline_open`·`terminal_lost` 분류
14. contact 입력 검증: 위·아래 전극 겹침, 같은 net, hole 하나가 전극 둘에 걸친 경우 오류
15. 자동 단자 open: 모서리 조각·단자 없는 가지 절단 무시, ROI 경계 단자, 명시 경로와의 동등성
16. CNO 완전·일부·무 가림, 적용 공정 외 조건, `contact_lost`
17. open·short·CNO 겹침에서 `1/k` 배분 합계, 전극·contact 위험도 합계

성능 검증은 기능 테스트와 분리해 동일 장비에서 warm-up 후 5회 측정한다. 중앙 TAT, 최소/최대,
변동계수, 최대 RSS, 출력 크기를 남긴다. 복합 200² failure mode별 측정과 같은 형상을 1000² 중앙에
둔 확장성 측정은 완료했다. 동일 형상의 공정 간 결과 재사용과 전극 bounding box 연결성 판정은
픽셀 동등성을 유지하면서 200² open을 5.27배 향상했다. padded 1000² 통합은 16.385초였지만 실제
제품의 긴 배선과 전극별 bounding box를 대표하지 않는다. 실제 제품형 입력이 확보되면 별도 측정한다.
산포는 무제한 open cache 없이 64조건을 먼저 측정한다.

## 6. 이번 보완에서 제외할 항목

다음은 코드 결함 수정과 기존 sim 대체 판정이 끝난 뒤 별도 사양으로 진행한다.

- TR의 실제 신호 분리 기능과 TR 손상 불량 판정
- 전도성/비전도성 이물별 접촉 거리와 전기적 short 조건
- 타원형, 길쭉한 형상, 불규칙 형상과 방향각 분포
- 실측 원본 이물 크기 표본의 경험분포 적분
- 실제 불량률/수율 변환과 Poisson 가정
- 동적 연결성 또는 bounded/disk cache를 사용하는 후속 open 고속화와 대형 실제 ROI 계산
- `pair_margin_um`의 선택 저장, 정밀도 또는 압축 정책 변경. 현재 계약과 수치 의미를 먼저 확정해야 한다.
- ~~CNO 가림 기준, 병렬 contact hole의 불량 기준~~ 2026-09-30까지 모두 확정(hole 중심 픽셀 덮음 · hole 하나라도 막히면 CNO)
- 같은 net의 다른 layer 우회 경로를 반영한 net 단위 open, 구동단·부하 정보에 따른 불량 심각도

다만 단계 2의 산포 mask 순차 처리는 실제 1000² 입력을 실행하기 위한 메모리 보완이므로 이번 P0에
포함한다. open 고속화는 정확도 동등성을 증명한 뒤 P1로 진행한다.

## 7. 최종 완료 조건

코드 보완 묶음은 다음 조건을 모두 만족해야 완료다.

- 전체 자동 테스트와 새 불변식 검증 통과
- 200²/1000² TAT·RSS 재측정 및 기존 기록과 비교
- schema 2.0과 예제 입력/출력 문서 일치
- clean commit 기준 재현 가능한 실행 manifest와 artifact hash 생성
- 기존 sim 실제 동일 입력 비교 및 미분류 차이 0
- 이후에만 실제 Fab 산포·평균 이물 자료와 simcompare 연결 진행

보고서와 PDF는 구현 중간마다 다시 만들지 않는다. 코드, 실제 비교 결과, 위험도 출력이 모두 고정된
뒤 마지막에 HTML과 PDF를 한 번 갱신한다.
