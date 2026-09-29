# 병합 준비 검토 — 2026-09-26

**최종 판정 (2026-09-27): 병합 가능 (PASS-WITH-NOTES).** 이 판정은 `fca4806`에 커밋하지 않은 수정 5개 파일을 더한 후보에만 적용된다. 수정본 diff 지문은 `13ae89c192d8b120d6d62fb39f760837509927c0a8600302c8f46fda6d5e0240`이다.
- F3은 사용자가 수락했다.
- F1·F2와 후속 V2-1·V2-2·V2-4는 수정했다.
- 새 Opus 5.5 high 독립 검토 v3의 판정은 PASS-WITH-NOTES다.
- 병합 전에 사용자가 수정본을 커밋해야 한다. 커밋 내용이 이 지문과 다르면 이 판정은 적용되지 않는다.

아래 1차 기록(v1 FAIL)은 이력으로 남긴다.

## 대상

| 항목 | 값 |
|---|---|
| 레포 | `/Users/jonghyeonlee/git/particle-sim` |
| 후보 브랜치 | `v5-extractor-integration` |
| 후보 HEAD | `fca4806fbd6eaa4f11139ea25e45b7d1188782ec` |
| 기준 main | `f4372519fa5aa61dc423ecd476f7224c845e6d6d` (= origin/main = merge-base, fast-forward 가능) |
| 커밋 수 | 19 (`main..HEAD`) |
| 제품 범위 | `src`, `tests`, `configs`, `requirements*.txt` 50개 파일, +7393 / −331 |
| 제품 범위 diff 지문 | sha256 `6c743e7fa789b9286a72b536985ee0380b9df9c15a53a00d5778164629e5cbb5` |

## Orca 기록

| 항목 | 값 | 결과 |
|---|---|---|
| Run | `run_84aae6435804` | 코디네이터 `term_ec27678b-…` (재개 세션, claude-opus-5-5) |
| Task 1 | `task_cd49739da79a` Astra 계획 (gpt-6-astra high, 화면 확인) | Dispatch `ctx_45a9d22e7a83` succeeded. 계획 v1 수락(수정 조건 포함) |
| Task 2 | `task_fe224e5a2992` 독립 검토 (claude-opus-5-5 high, 새 세션, 화면 확인) | Dispatch `ctx_faad8f19d8c1` succeeded. 판정 FAIL |
| 워커 터미널 | `term_7cec63e1-…`(Astra), `term_22bfb28f-…`(검토자) | release 결과 `retained`/`external_terminal`. 명시 터미널이어서 Orca가 닫지 않는다. reclaimable 0 |

증거 파일은 모두 git-ignored `outputs/merge-review-20260926/`에 있다.

- 계획과 수락: `astra-plan-v1.md`, `plan-v1-acceptance.md`
- 코디네이터 검증: `validation-v1/`
- 독립 검토 보고서와 프로브: `review-v1/independent-review.md`, `review-v1/*.py|json`

## 코디네이터 검증 (후보 fca4806)

| 검증 | 결과 |
|---|---|
| pytest 전체 | **369 passed** (51 s) |
| ruff `check --no-cache src tests` | 통과 |
| 의존성 (`uv pip check`) | 43개 패키지 모두 호환 |
| main 테스트를 HEAD `src`로 실행 | 253 passed, 1 failed. `test_union_mismatch`는 의도된 완화다(F6 참고) |
| V2 main 대 HEAD (schema-1 합성 4종과 probability) | 삭제 0건. 공통 키 값 차이 0건. 출력은 추가만 있다. `schema_version` 1.0→2.0. 같은 seed 재실행 결과 동일. 이전 `apply-fab --ca` 결과 동일. 음수 경로는 exit 2이고 출력이 없다 |
| V3 contact/CNO | 변동 0 probability가 nominal 지도와 같다. 작은 ROI 분석은 exit 2이고 출력이 없다. 이전 기록의 "PB 0 px"는 step 축 합계였다. unit별로는 `[412, 64, 64]`이므로 결함이 아니다(F13) |
| F2 재현 (코디네이터) | main이 만든 schema-1 probability run에 HEAD `make-report`를 실행하면 exit 1(`IndexError`)이고, main은 exit 0이다 (`validation-v1/f2-repro/`) |

## 독립 검토 결과 (Q1–Q7)

| Q | 판정 | 요지 |
|---|---|---|
| Q1 확률·CA·배분·Fab·보존 | PASS | 1/k 배분을 손계산과 대조해 정확히 일치했다. 전극·layer 합이 보존되고 Fab 밀도 의미가 정확하다 |
| Q2 fast open = 전체 판정 | PASS | 모든 ROI 중심에서 독립 oracle과 비교해 불일치 0건. 단조성 반례도 유지된다 |
| Q3 CNO 규칙 일관성 | **FAIL** | 규칙 자체는 oracle과 0건 불일치다. unit ID 충돌은 F1 |
| Q4 schema-1 호환 | **FAIL** | F2, F3 |
| Q5 결정성·packed | PASS | memory와 packed 경로 결과가 비트 단위로 같다. seed 중복 제거, 임시 디렉터리 정리 |
| Q6 오류 경로·무결성 | PASS | Minor F4, F5, F10, F11 있음 |
| Q7 게이트·잔여 위험 | **FAIL** | 게이트는 통과했지만 Major가 남아 있다 |

## 차단 항목 (병합 전 필수)

- **F1 Major — CNO unit ID 충돌.** `geometry/contacts.py:cno_unit_id`가 `f"cno:{contact}:{lower}__{upper}"`로 ID를 만든다. 전극 ID에 `__`가 있으면 서로 다른 링크가 같은 ID가 된다(`P__Q`+`W` = `P`+`Q__W`). 이때 `build_contact_topology`의 `links.setdefault`가 두 링크를 조용히 합친다. 결과는 CNO topology와 전극 귀속 오류다. 수정 방향은 두 가지다. contacts가 선언되면 전극 ID의 `__`를 거부하거나, 같은 unit_id에 다른 (lower, upper)가 오면 오류를 낸다. 회귀 테스트도 추가한다.
- **F2 Major — `make-report` 회귀.** `visualization/hotspots.py:plot_convergence` 159행은 `scramble_seed` 열이 없는 schema-1 convergence에서 1-튜플 키에 `key[1]`을 적용해 `IndexError`를 낸다. `reports.py:make_probability_report`도 v1에 없는 `ca_union_*` 열을 읽는다. 수정 방향은 두 가지다. v1 run을 지원하거나, v1 run을 `InputContractError`(exit 2)로 명확히 거부한다.
- **F3 Major (결정 항목) — 출력 schema 2.0으로 무조건 전환.** `io/outputs.py:SCHEMA_VERSION = "2.0"`. V2 결과로는 추가만 있는 변경이다. 하지만 `schema_version == "1.0"`을 검사하는 외부 소비자는 새 run을 거부한다. 제품 책임자가 이 전환을 명시적으로 수락하고, 지원하지 않는 v1 소비자를 기록해야 한다. schema-1 입력에 대한 변환 경고 여부도 함께 결정해야 한다.

## 비차단 항목 (후속 권장)

| 항목 | 등급 | 내용 |
|---|---|---|
| F4 | Minor | `analyze-central --run <probability run>`은 config_hash 불일치 때문에 항상 실패한다. 이 분기는 실행되지 않는 코드다 |
| F5 | Minor | `apply-fab --run`은 `product_id`만 확인하고 design_version이나 manifest hash는 확인하지 않는다 |
| F6 | Minor | 전극 ⊆ layer 완화가 schema-1에도 적용된다. 할당되지 않은 layer 픽셀이 있어도 경고가 없다 |
| F7 | Minor | 기존 문제. pair_id와 path_id가 같아도 허용되어 결과가 합쳐진다 |
| F10 | Minor | legacy 동일성 확인이 mask hash 일부만 맞아도 통과한다 |
| F8, F9, F11–F16 | Note | 실패 run의 부분 산출물(main과 같음), git provenance가 cwd를 따른다, dedup이 raw digest를 확인하지 않는다, `auto_open.py:229` EOF 공백, 의도된 schema-1 강화 2건, 주석 부정확, scenario 제외 의미 |

## 2차: 사용자 결정과 F1·F2 수정 (2026-09-26~27)

- **F3은 수락한다.** 출력 schema를 무조건 2.0으로 바꾸는 것을 사용자가 2026-09-26에 받아들였다("f3는 받아들이고 나머지 수정 진행해"). 코드는 바꾸지 않는다.
- **수정 Task**: `task_d080f51c1120`. Dispatch는 `ctx_5eda2b735a69`, 터미널은 `term_908346e6-…`(gpt-6-sol medium, 화면 확인, Codex 0.157.1)이다. 결과는 succeeded이고 코디네이터가 수락했다.
  - F1: `build_contact_topology`에서 같은 unit_id에 다른 (contact, lower, upper)가 오면 `InputContractError`를 낸다. strict와 non-strict 모두 적용되고 CLI exit는 2다. 구분자 모호성을 일반적으로 막으면서, 이미 쓰이는 전극 ID는 그대로 허용한다.
  - F2: `plot_convergence`에 1-튜플 가드(`len(key) > 1`)를 넣었다. `make_probability_report`는 `pair_ids`로 대체하고, v1에 없는 union 통계는 "해당 없음"으로 표시한다.
- **수정 후보**: fca4806 + 커밋하지 않은 수정 5개 파일(+84/−5). `git diff fca4806 | shasum -a 256` = `da9679ae2c7d8729a30fada3a736e7acf4af0f9d7bbab381629d1eee0bbb40ce`
- **코디네이터 재검증** (`outputs/merge-review-20260926/validation-v2/`)
  - pytest 전체 **373 passed**. ruff 통과, `git diff --check` 통과.
  - F2 재현: main이 만든 schema-1 probability run에 대해 `make-report`가 exit 0으로 끝나고, 보고서에 "해당 없음"이 표시된다.
  - 새 회귀 테스트를 fca4806 원본 src에서 실행했다. F1은 strict·non-strict 모두 `DID NOT RAISE`로 실패하고, F2 legacy는 `IndexError`로 실패한다. schema-2 케이스는 통과한다. 수정본에서는 모두 통과한다.
  - schema-2 보고서가 바뀌지 않았다. 같은 HEAD run에 대해 base와 수정본의 `report.md` diff가 0줄이고 PNG 5개가 바이트 단위로 같다.
- **독립 검토 v2**: `task_ccee49402f7c`, Dispatch `ctx_056913f25d06`, 터미널 `term_ccb47cfd-…`(claude-opus-5-5 high, 새 세션, 화면 확인). 판정 **FAIL, Major 1건**(V2-1)이다.
  - V2-1: `units.py:evaluate_units`는 변형된 topology를 nominal `unit_id`로 조회한다. 그래서 GI overlay가 hole을 충돌하는 `(P, Q__W)` 쌍에 다시 연결하면, 그 결과가 nominal unit `(P__Q, W)`에 조용히 귀속된다. 코디네이터가 코드와 프로브로 확인했다.
  - V2-2 Minor: v1 union 평균이 nominal 열에 표시된다.
  - V2-3 Minor: 실행 도중 중단된다.
  - V2-4 Note: 오류 메시지가 영어다.
  - V2-5 Note: 테스트가 변형 경로를 다루지 않는다.

## 3차: V2-1 후속 수정과 최종 검토 (2026-09-27)

- **후속 수정 Task**: `task_b187b17ae0c4`, Dispatch `ctx_115f6c54825d`. 같은 Sol 터미널(`term_908346e6-…`)을 재사용했다. 결과는 succeeded이고 코디네이터가 수락했다.
  - V2-1: strict 로드 시 contact마다 `resolve_contact_layers`로 정한 아래·위 layer에서 net이 같은 모든 전극쌍을 열거한다. 서로 다른 쌍이 같은 `cno_unit_id`를 만들면 `InputContractError`(exit 2)를 낸다. 코디네이터는 이 열거가 non-strict 경로에서 만들어질 수 있는 링크 집합과 정확히 같다는 점을 코드로 확인했다. 두 경로는 layer 기준이 같고, 다른 net은 `lost` 처리된다. 또 `contact_id`에는 `:`와 `__`가 금지되어 있다. 이 변경으로 V2-3도 없어진다. 기존의 topology 내부 검사는 이중 방어로 유지한다.
  - V2-2: v1 행에서 v1 평균을 `ca_union_mean_um2`에, `ca_pair_sum_nominal_um2`를 `ca_raw_unit_sum_mean_um2`에 표시한다. nominal, p05, p95는 "해당 없음"이다.
  - V2-4: 충돌 메시지를 한국어로 바꿨다.
- **최종 후보**: fca4806 + 커밋하지 않은 5개 파일(+167/−5). 지문은 `13ae89c192d8b120d6d62fb39f760837509927c0a8600302c8f46fda6d5e0240`
- **코디네이터 재검증** (`outputs/merge-review-20260926/validation-v3/`)
  - pytest 전체 **375 passed**. ruff 통과, `git diff --check` 통과.
  - main이 만든 schema-1 probability run에서 `make-report`가 exit 0으로 끝나고, v1 평균이 mean 열에 표시된다.
  - schema-2 보고서는 fca4806 출력과 같다. `report.md` diff가 0줄이고 PNG 5개가 바이트 단위로 같다.
  - 새 테스트 4개(F1 strict·non-strict, 시나리오에서만 생기는 충돌, F2 legacy)는 fca4806 원본 src에서 실패한다.
- **독립 검토 v3**: `task_42dee4397964`, Dispatch `ctx_7e3ec17a4f4b`, 터미널 `term_9cd95261-…`(claude-opus-5-5 high, 새 세션, 화면 확인). 판정 **PASS-WITH-NOTES**. 보고서는 `outputs/merge-review-20260926/review-v3/independent-review.md`다.
  - S1~S4 모두 PASS, S5는 PASS-WITH-NOTES다.
  - V2-1 geometry는 `validate-input`·`simulate-nominal`·`simulate-probability`에서 로드 시 exit 2로 끝나고 출력을 남기지 않는다.
  - 정상적인 `__` ID, 다른 net, layer를 바꾼 manifest는 모두 로드된다. overlay를 여러 값으로 바꿔 봐도 spec과 link가 어긋나는 경우가 없다.
  - 새 테스트는 fca4806과 이전 수정본에서 실패하고, 최종 후보에서 40/40 통과한다.
  - 남은 것은 Note 2건뿐이다. V3-1: contact마다 O(L·U) 열거를 한다(같은 net 전극이 layer당 1000개일 때 0.27 s). V3-2: 정적 거부는 의도적으로 보수적이다.

## 최종 판정 근거와 후속 항목

- 병합 판정: **병합 가능 (PASS-WITH-NOTES)**. Major 3건 가운데 F1·F2는 수정했고 F3은 사용자가 수락했다. v2에서 나온 Major도 닫혔다. Blocker, Major, Minor는 남지 않았고 v3 Note 2건만 있다.
- 병합 조건
  1. 사용자가 수정본 5개 파일을 커밋한다. 커밋 diff가 지문 `13ae89c1…`과 같아야 한다.
  2. fast-forward 병합한다. main이 `f437251`에서 움직이면 다시 확인한다.
- 후속 항목(이번 범위 밖): v1의 F4, F5, F6, F7, F8, F10(Minor)과 F9, F11–F16(Note), v3의 V3-1과 V3-2(Note).
- **운영 기록**: 첫 Sol 터미널에서 Codex 업데이트 프롬프트가 떴다. 코디네이터가 Skip을 보내려 했지만 입력이 기본값(Update now)으로 들어가, Codex가 0.157.0에서 0.157.1로 전역 업데이트되었다. 이 터미널에는 Dispatch가 없었으므로 닫고 같은 명령으로 다시 띄웠다. `~/.codex/hooks.json`은 바뀌지 않았음을 확인했다.

## 보존 상태

- 제품 파일 수정, 커밋, 병합, push, 브랜치 변경은 하지 않았다.
- 추적되지 않은 `docs/layout-electrode-extraction-plan.md`는 보존했다.
- 이 문서는 당시 커밋하지 않은 새 파일이었다(2026-09-29 커밋).

## 커밋과 병합 (2026-09-27, 사용자 지시)

- 수정 5개 파일을 `3290d3e`로 커밋했다. `git diff fca4806 3290d3e -- src tests`의 sha256은 `13ae89c1…0240`이다. 검토한 후보와 같으므로 PASS-WITH-NOTES 판정이 이 커밋에 적용된다.
- 원격과 로컬 main이 `f437251`인 것을 확인한 뒤, 로컬 main을 `3290d3e`로 fast-forward했다. 병합 커밋은 없다.
- push는 하지 않았다. 원격 main은 아직 `f437251`이다.

## Minor 후속 수정 (2026-09-29)

브랜치 `fix/merge-review-minors`(기준 main `29a21e0`)에서 v1 Minor 5건을 고쳤다. 각 회귀 테스트는 수정 전 `src`에서 실패하고 수정본에서 통과한다. 전체 pytest 821 passed, `ruff check src tests` 통과.

| 항목 | 수정 | 회귀 테스트 |
|---|---|---|
| F4 | `analyze-central --run`이 probability run이면 `config_hash`(scenario 포함) 대신 `design_version`, `input_hashes.manifest`, mask hash로 같은 제품인지 확인한다. nominal run은 기존 `config_hash` 대조를 유지한다 | `test_central_window.py::test_probability_run_is_matched_by_manifest_hash` |
| F5 | `apply-fab --run`이 `design_version`과 제품 YAML의 sha256(`input_hashes.manifest`)도 대조한다. 다르면 출력 폴더를 만들기 전에 exit 2 | `test_cli_phase2.py::test_apply_fab_run_rejects_other_design_or_manifest` |
| F6 | 전극에 할당되지 않은 layer 픽셀이 있으면 schema와 무관하게 `run_manifest.warnings`에 `LAYER_UNASSIGNED:` 경고를 남긴다(nominal·probability 모두). 판정 동작은 바꾸지 않았다 | `test_io_masks_stack.py::test_layer_may_contain_unassigned_tr_region`, 합성 실행에서 경고 기록 확인 |
| F7 | manifest 검증에서 `pair_id`와 `path_id`가 겹치면 거부한다 | `test_schema_config.py::test_pair_and_open_path_ids_must_differ` |
| F10 | `compare-sim` 기준 YAML의 `input_hashes`가 현재 실행의 모든 mask hash(`manifest`·`scenario` 제외)를 포함해야 한다. 예시 YAML과 입출력 안내서를 함께 고쳤다 | `test_legacy_sim_compare.py::test_reference_from_another_product_or_mask_is_rejected` |

호환성 영향: F5는 실행 뒤 제품 YAML을 수정하면(예: `fab_id` 변경) `apply-fab --run`을 거부한다. 이때는 수정한 YAML로 다시 실행해야 한다. F10은 일부 mask만 적은 기존 비교 YAML을 거부한다.

남은 항목: F8, F9, F11–F16(Note), V3-1, V3-2.

독립 검토(새 Opus 세션, high, `5dccfc6` 대상): **PASS-WITH-NOTES**. Minor 1건(F4 테스트가 mask hash 불일치를 검사하지 않음)과 F6 호출부 자동 테스트 부재(Note)는 테스트를 추가해 보완했다. 나머지 Note는 다음과 같다.
- F6 개수는 scenario 변형 전 명목 mask 기준이다.
- F10은 terminal·contact mask hash까지 요구한다. 기존 sim이 읽지 않은 파일이어도 run_manifest에서 복사해 적어야 한다.
- F5 때문에 실행 뒤 제품 YAML을 수정하는 흐름은 막힌다(의도된 동작).

같은 정리에서 추적되지 않던 `docs/layout-electrode-extraction-plan.md`를 삭제했다(커밋 diff에는 나타나지 않는다). `layout-electrode-extractor` 저장소의 `docs/PLAN.md`와 바이트 단위로 같고 그 저장소에서 관리한다.
