# 보고서 검토와 후속 계획

검토일: 2026-09-18. 대상 HEAD: `f4372519fa5aa61dc423ecd476f7224c845e6d6d`.
범위: 보고서 전체 17절, 인쇄 결과, 저장된 실행 근거, 관련 구현·계약 대조.
전체 코드 감사나 실제 제품의 물리적 유효성 검증은 아니다. 아래 내용 대조는 당시 PDF 근거를
기록한 것이며, 2026-09-19에 P0 제품 코드 보완을 별도로 완료했다.

구현 변경 항목의 현재 동작, 수정 이유, 대상 파일과 완료 기준은
[`../CODE-CHANGE-BACKLOG.md`](../CODE-CHANGE-BACKLOG.md)에 모아 두었다.

## 2026-09-19 구현 갱신

- layer mask 부분집합·잔여 TR 통계, 전극 조합 감사, 단자 영역 open과 적용 공정을 구현했다.
- x/y overlay와 항목별 근거, 명목/산포 union CA 분리, 여러 Sobol seed와 packbits 저장을 구현했다.
- Fab 평균 이물 개수와 검사면적 입력, layer·전극 조합·open·개별 전극 위험도와 manifest를 구현했다.
- 외부 기존 sim 비교 엔진과 예제 계약을 구현했고 전체 **279 tests**가 통과했다.
- 복합 open 결과의 공정 간 재사용과 전극 bounding box 연결성 판정을 적용해 200² open을 최초
  기준선 대비 5.27배 개선했고, hotspot·CA 동등성을 확인했다.
- 같은 국부 형상을 padded 1000²에서 측정한 통합 TAT는 16.385초, 최대 RSS는 843.6MB,
  출력은 621.4MB였다. 이는 실제 긴 배선 대표값이 아니라 ROI·출력 확장성 값이다.
- 실제 현장 sim 산출물이 없어 동일 입력 XOR/CA 대조와 대체 판정은 남아 있다.

## 2026-09-20 구현 갱신

- 미커밋 작업 트리 검토에서 찾은 결함 8개를 모두 수정했다. open은 이물 없는 상태의 연결을 먼저
  확인하고, 이미 끊긴 경로와 사라진 단자를 이물 기인 CA에서 제외한다.
- contact·via mask 입력을 추가했다. contact hole과 ROI 경계를 open 자동 단자로 쓰고, contact가 막히는
  CNO(Contact Not Open)를 별도 불량 유형으로 계산한다. CNO 가림 기준은 잠정(전부 가림)이다.
- 기존 sim 비교 계약을 schema 2.0으로 올려 제품·설계 버전·mask hash가 같아야 비교한다.
- 복합 200² 파일럿에서 기존 결과는 픽셀 단위로 유지됐고, 단자를 contact로 바꾼 자동 단자 open은
  명시 단자 결과와 같았다(막힌 contact 위치만 CNO로 이동). 전체 테스트는 303개가 통과했다.
- 보고서에 추가할 항목: CNO 불량 유형과 잠정 기준, 자동 단자 open, 입력 계약 변경, 아래 TR 개선 항목.
  사용자 요청으로 2026-09-20에 HTML/PDF를 25쪽으로 다시 만들었다. 불량 유형별 장(4~10쪽)과
  생각해야 할 부분(23~25쪽)을 추가했고, 넘침·글자 손실 검사를 통과했다(`evidence/pdf-validation.json`).
  `REPORT-SCRIPT.md` 발표 원고는 새 장 구성을 아직 반영하지 않았다.

## 2026-09-21 발표 원고 갱신

- `REPORT-SCRIPT.md`를 25쪽 판 기준으로 다시 썼다. 쪽마다 핵심 메시지·발표 스크립트·사실과 한계·첨삭
  요청을 두고, 판정이 갈릴 수 있는 9개 쪽에는 예상 질문을 덧붙였다. 17쪽 판과의 대응표와 발표 시간
  배분도 넣었다.
- 원고의 모든 수치는 `report.html` 본문과 대조해 일치를 확인했다. 인용한 코드 식별자
  (`CNO_COVERAGE_RULE`, `cno_coverage_rule`, `risk_by_contact.parquet`, `insufficient_terminals`,
  `baseline_open`, `terminal_lost`, `contact_lost`, `new_unit`)도 현재 소스에 존재한다.
- 이전 17쪽 판 원고는 `git show 9013b9d:docs/report/REPORT-SCRIPT.md`에 있다.
- 첨삭은 원고가 아니라 `build_report.py`에 반영한 뒤 HTML/PDF를 한 번 다시 만든다. 이번 갱신에서는
  보고서 본문·그림·PDF를 바꾸지 않았다.

## 이번에 확인·수정한 사항

- 기존 PDF는 17쪽이며, 9·10·11·16·17쪽에 앞 절의 내용이 넘어왔다. A4 인쇄 폭에서 본문이 고정 높이를 초과하고 절대 위치 footer와 충돌했다. 표의 분할 방지 규칙까지 작용해 다음 절과 섞였다.
- 인쇄 본문 폭을 178 mm로 고정하고, 표·그림·문단 간격과 긴 부록의 코드 행간을 조정했다. 본문을 숨기거나 삭제하지 않았다.
- Chrome에서 폰트·이미지를 기다린 뒤 출력하도록 변경했다. 출력 전 footer 여유·수평 넘침·이미지 로드, 출력 후 페이지 수·쪽별 텍스트 문자 수량·PDF 텍스트 경계를 검사한다. 실패하면 기존 PDF를 교체하지 않는다.
- 재생성 결과는 `evidence/pdf-validation.json`에 기록하고, 페이지 이미지는 `work/page-01.png` 등에 보존한다.
- 기존 PDF는 `work/particle-sim-report-before.pdf`에 보존했다.

## 내용 대조 결과

| 보고 내용 | 대조 근거 | 판단 |
|---|---|---|
| pilot d=10의 pair CA 및 union 1750 / 500 / 0 µm² | evidence/run/ca.parquet, ca_union.parquet | 일치 |
| 확률 CA 650.78125, p05 557.5, p95 750 µm² | evidence/prob/ca.parquet | 일치 |
| common 0.0004416015625 / actual 0.00054459375 | evidence/prob/fab_risk.parquet의 모드별 expected_defect_count 합 | 일치, pair 합이며 union 위험 아님 |
| triple_overlap 10 µm: union 350, pair 합 650 | evidence/triple/ca_union.parquet | 일치 |
| 257 passed / 44.91초 | evidence/pytest.log | 복합 fixture 추가 후 전체 테스트 실행 기록과 일치 |
| 279 passed / 44.76초 | evidence/pytest-p0.log | P0 보완과 open bbox 회귀 테스트를 포함한 현재 전체 실행 |
| 산포 0·seed 재현·open 전수 대조 | evidence/recheck.json 및 recheck.py | 저장된 검증 범위와 일치 |
| 기존 방식 차이 3 / 30 / 20 px | evidence/regression-probes.json | 일치 |
| 당시 open 최적화 미반영 | 이전 failure/open_failure.py, units.py | 당시 전체 ROI label 사용. 2026-09-19 공정 간 cache와 전극 bbox 판정으로 보완 |
| simcompare의 vertical Fab 공정 충돌 | export/simcompare.py 및 인접 저장소 data_contract.md | 계약 불일치가 남음 |

원 사용자 원자료 전체와의 요구사항 추적표는 없으며, 이번 검토는 저장소 계약 및 보존된 근거 중심이다. 원자료 PDF 추출문은 work/source-pdfs.txt에 있지만 원본 일부와 하위 계약은 외부 경로를 참조한다. 보고서 폴더만 전달하면 모든 출처를 독립적으로 재검토할 수는 없다.

## 추가 작업 계획

추가 작업이 있으므로 PDF 재생성만으로 제품 개발 완료로 판단하지 않는다. P0 구현은 완료했으며
아래 표는 현재 상태로 갱신했다.

| 순서 / 우선순위 | 작업·담당 | 산출물 및 완료 기준 | 선행 조건 |
|---|---|---|---|
| 완료 | P0 코드 보완·비교 엔진 | schema 2.0, 위험도 출력, provenance, 메모리 보완, compare-sim과 279 tests | 실제 비교용 코드는 준비됨 |
| 완료 | 검토 결함 수정·contact/CNO · 개발 | open 기준 연결 검사, 자동 단자 open, CNO, 비교 계약 2.0, 303 tests, contact 파일럿 동등성·TAT | CNO 가림 기준과 병렬 hole 기준은 논의 필요 |
| 다음 / P0 | 기존 `sim` 비교 입력 준비 · 개발/공정 | 기존 `sim`의 입력·출력, 좌표축, 픽셀 크기, 직경, 이물 유입 공정과 전극 조합 대응표 확정 | 기존 `sim` 실행 환경과 대표 입력·출력 확보 |
| 다음 / P0 | 기존 `sim` 대조와 대체 판정 · 개발/공정 | 공통 전극 조합의 hotspot 겹침·XOR·CA 차이를 비교하고 모든 차이를 분류. 기존 계산 대상이 아닌 G–A 같은 대각 조합은 `N/A/new_pair`로 보고. 미분류 차이가 없고 대표 ROI의 TAT/RSS가 허용 범위일 때 대체 후보로 판정 | 동일한 mask·축척·공정·크기 조건 |
| 4 / P1 | 실제 제품 기준·산포 pilot · 개발/공정 | 실제 mask에서 기준 short를 먼저 확인한 뒤 mask layer별 Fab 실측 overlay·CD만 독립 입력으로 적용. 같은 seed의 Sobol 64→128→256조건과 여러 scramble seed로 전극 조합별·중복 제거 전체 CA 수렴과 TAT 확인 | 기존 `sim` 대체 판정. overlay 기준 layer·x/y·CD 정의와 원자료 ID 확정 |
| 5 / P1 | 이물 Fab 분포의 사후 적용 · 공정/개발 | 크기별 시뮬레이션을 완료한 뒤 공통 이물 분포를 union CA에 적용해 구조 상대 비교. 실제 Fab 공정·크기 구간 연결이 확정되면 실제 분포 결과를 별도로 계산 | 이물 유입 공정↔Fab 측정 공정, 대표 크기↔Fab bin 계약 |
| 6 / P1 | `sim-compare` 연결 · 양 도구 담당 | 이물 유입 공정, Fab 측정 공정과 상·하부 전극 layer를 분리한 계약으로 내보내기→읽기→가중 결과를 수기 계산과 대조 | `particle-sim` 대체 결과와 union/Fab 계약 확정 |
| 완료 | open 기준선과 1차 최적화 · 개발 | 복합 200² mode별 5회 측정, 공정 간 결과 재사용과 전극 bbox 판정 적용. open 5.27배·통합 4.04배 개선 및 픽셀/CA 동등성 확인 | 실제 제품형 긴 배선의 별도 측정은 남음 |
| 8 / P1 | 폐쇄망 반입 재현 · 개발/인프라 | 대상 OS·CPU·Python용 wheelhouse로 빈 환경에서 오프라인 설치, 전체 시험과 대표 계산 실행. 설치 로그·시간·hash·버전 보존 | 대상 플랫폼과 uv 또는 pip 절차 확정 |
| 9 / P2 | 대형 ROI와 추가 물리 모델 · 개발/공정 | 실제 최대 ROI에서 메모리 한계를 확인한 뒤 tile+halo 검토. TR 전기 기능, 이물 재질, 비원형 형상과 예상 불량 변환은 각각 물리 근거와 완료 기준을 마련한 후 별도 개발 | 실제 성능 병목과 측정·물리 자료 확보 |

현재 최우선 순서는 **기존 `sim` 대표 입력·출력 확보→동일 입력 비교→대체 판정**이다.
`sim-compare` 연결은 기존 `sim` 대체 판정과 결과 계약이 확정된 뒤
진행한다. 실제 제품 산포에는 측정된 mask별 overlay·CD만 사용하고, 근거가 없는 산포는 0으로
고정한다. 이물 크기별 Fab 분포는 접촉 시뮬레이션 후 적용한다.

## PDF 제작 시 추가할 향후 개선 항목

아래 내용은 현재 구현 범위가 아니며, 다음 PDF 재생성 때 "향후 개선" 절 또는 9번 작업(추가 물리
모델)의 설명으로 추가한다. 현재 결과의 의미를 바꾸지 않도록 구현 전 물리·출력 계약을 먼저 확정한다.

### TR에 대한 이물 영향 (2026-09-19 정리)

현재 TR은 layer mask에서 전극에 할당되지 않은 영역으로 픽셀 수·면적만 기록하며, short·open·hotspot
판정에는 포함하지 않는다. 반도체 분야의 Inductive Fault Analysis(IFA)는 이물 통계로 물리 결함을
만들고, 이를 트랜지스터 수준 결함(bridge, cut, gate-oxide short)으로 추출한 뒤 결함 유형별 critical
area로 발생 가능성을 순위화한다. 이 흐름을 현재 엔진에 맞춰 세 단계로 확장하는 방안을 제안한다.

**1단계 · TR 피격을 기하 위험으로만 기록 (우선 권장)**

- manifest에 `transistors` 항목을 추가한다. TR ID, 채널(active) 영역 mask, gate·source·drain 전극 ID를
  입력한다.
- 이물 원판이 채널 영역과 겹치는 중심 위치를 `tr_hit` hotspot·CA로 별도 출력한다. open 후보 계산과
  같은 원판 겹침 판정을 사용한다.
- 물리 의미를 확정하기 전에는 `union:any`에 합치지 않고 별도로 보고한다.

**2단계 · TR을 스위치로 보고 결함 유형 규칙으로 판정**

TR의 신호 분리 기능은 source와 drain이 다른 net이며 채널을 통해서만 조건부로 연결된다는 뜻이다.
이 관점에서 결함을 다음처럼 나누면 기존 엔진을 재사용할 수 있다.

| TR 결함 | 기하 조건 | 재사용 엔진 | 전기적 의미 |
|---|---|---|---|
| S–D short (항상 켜짐) | 도전성 이물이 채널을 가로질러 S·D에 동시에 접촉 | lateral short (S–D 조합) | 신호 분리 기능 상실 |
| gate short | 이물이 게이트 절연막을 관통해 G와 S/D 또는 채널에 접촉 | vertical short (G–S, G–D, 높이 차 포함) | 게이트 누설 또는 short |
| 채널 open (항상 꺼짐) | 채널 영역이 S 접점과 D 접점 사이에서 끊김 | open (채널 mask를 전극, S·D 접점을 단자 영역으로) | 구동 불가 |

각 규칙에는 open 경로와 같은 적용 이물 유입 공정과 이물 재질(도전성·비도전성) 조건이 필요하다.
재질별 접촉 판정 후속 항목과 같은 입력 계약으로 묶는다.

**3단계 · 회로·화소 불량 연결과 kill ratio 보정**

- TR 결함이 화소·라인 불량으로 이어지는 관계는 화소 회로와 netlist에 따라 다르다. 채널 피격이
  Vth 이동에 그칠 수도 있다.
- 결함 유형별 kill ratio(피격이 실제 불량이 되는 비율, 1 이하)를 두고 검사 좌표와 전기 검사·FA
  결과를 대조해 보정한다.
- 보정 전까지 출력은 "기하 피격 기대값"으로 유지하고 예상 불량률로 표기하지 않는다.

**구현 전 확정할 사항**

1. TR 종류와 구조(예: LTPS/Oxide TFT, top/bottom gate). vertical 판정의 높이 순서를 정한다.
2. 채널(active) 영역 mask와 S/D 접점 위치를 입력으로 제공할 수 있는지
3. `tr_hit`을 `union:any`와 분리 보고할지
4. kill ratio 보정에 쓸 불량 좌표·FA 자료 확보 여부

참고 문헌:

- [Carafe: An Inductive Fault Analysis Tool for CMOS VLSI Circuits](https://www.researchgate.net/publication/2705927_Carafe_An_Inductive_Fault_Analysis_Tool_for_CMOS_VLSI_Circuits)
- [Automatic Fault Extraction and Simulation of Layout Realistic Faults](https://www.cs.york.ac.uk/rts/docs/SIGDA-Compendium-1994-2004/papers/1995/edt95/pdffiles/09c_2.pdf)
- [Systematic Extraction of Critical Areas From IC Layouts](https://link.springer.com/chapter/10.1007/978-1-4757-9957-6_4)
- [IC spot-defect and fault semantics — a unified framework](https://www.researchgate.net/publication/224734988_IC_spot-defect_and_fault_semantics-a_unified_framework)
- [Modeling Gate Oxide Short Defects in CMOS Minimum Transistors](https://www.researchgate.net/publication/3965601_Modeling_Gate_Oxide_Short_Defects_in_CMOS_Minimum_Transistors)
- [System and method for classifying defects in an electrical circuit (TFT array 결함 분류)](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6982556)
- [TFT Backplane and Issues for OLED](https://www.researchgate.net/publication/326869223_TFT_Backplane_and_Issues_for_OLED)

이 조사는 방법론 수준이며, 현재 제품 공정에 맞는 TR 판정 규칙은 공정 담당 확인이 필요하다.

## PDF 재생성

저장소 루트에서 실행한다. 계산용 Python 환경과 PDF 도구 환경을 분리한다.

```sh
.venv/bin/python docs/report/build_report.py
uv venv docs/report/work/pdf-tools
uv pip install --python docs/report/work/pdf-tools/bin/python playwright==1.63.0 pymupdf==1.28.2
zsh docs/report/render-pdf.sh
```

이미 PDF 도구 환경이 있으면 환경 생성·설치를 생략한다. 로컬 Google Chrome을 사용하며 브라우저를 따로 내려받지 않는다. 다른 환경은 `REPORT_PYTHON`에 도구가 설치된 Python 경로를 지정한다. 검증용 이미지는 matplotlib 그림 원본의 라벨 의미까지 자동 검증하지 않으므로 최종 육안 확인을 병행한다.
