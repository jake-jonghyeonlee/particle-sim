from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parent
E = R / "evidence"
pages: list[str] = []


def table(headers, rows, cls=""):
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def fig(name, caption, cls=""):
    return (
        f'<figure class="{cls}"><img src="figures/{name}" alt="{html.escape(caption)}">'
        f"<figcaption>{caption} <b>합성 데이터.</b></figcaption></figure>"
    )


def note(text):
    return f'<div class="note">{text}</div>'


def source(text):
    return f'<div class="source">근거 · {text}</div>'


def page(kicker, title, body):
    number = len(pages) + 1
    pages.append(
        '<section class="page">'
        '<header>PARTICLE-SIM <span>기술개발 내부자료 · 개발 결과 보고서</span></header>'
        f'<div class="kicker">{kicker}</div><h1>{title}</h1>{body}'
        '<footer><span>2026.09.19 · 기하 접촉 위험 / 합성 검증</span>'
        f'<b>{number:02d}</b></footer></section>'
    )


large_short = json.loads((E / "large-short.json").read_text(encoding="utf-8"))
large_prob = json.loads((E / "large-prob.json").read_text(encoding="utf-8"))
large_prob_p0 = json.loads((E / "large-prob-p0.json").read_text(encoding="utf-8"))
failure_mode_pilot = json.loads(
    (E / "benchmark-failure-modes-200.json").read_text(encoding="utf-8")
)
mode_tat = failure_mode_pilot["tat_summary"]
failure_mode_pilot_final = json.loads(
    (E / "benchmark-failure-modes-200-after-open-crop.json").read_text(encoding="utf-8")
)
mode_tat_final = failure_mode_pilot_final["tat_summary"]
padded_1000 = json.loads(
    (E / "benchmark-failure-modes-padded-1000.json").read_text(encoding="utf-8")
)["tat_summary"]
mode_ca = {
    row["report_mode"]: row
    for row in failure_mode_pilot["combined_result"]["mode_contribution"]
}
complex_validation = json.loads((E / "complex-validation.json").read_text(encoding="utf-8"))
bench = {
    name: json.loads((E / f"benchmark-{name}.json").read_text(encoding="utf-8"))["summary"]
    for name in ("pilot-200", "complex-200", "pilot-1000", "complex-1000")
}
triple = pd.read_parquet(E / "triple/ca_union.parquet")

page(
    "01 / 개발내용 요약",
    "다층 전극의 이물 위험 계산을 구현하고<br>합성 형상으로 검증했다",
    '<p class="lead">particle-sim은 이물 중심 후보 위치마다 <b>lateral short, vertical short, open</b> 가능성을 계산하고, 그 위치의 면적을 적분해 Critical Area(CA)를 구하는 2.5D 시뮬레이터다.</p>'
    + note('<b>현재 단계 · P0 코드 보완 완료</b><br>입력·산포·Fab·출력 계약과 외부 기존 sim 비교기를 구현했다. 실제 제품 mask와 기존 현장 sim 산출물의 동일 입력 대조가 다음 단계다.')
    + '<div class="metrics"><div><strong>279 / 279</strong><span>전체 pytest 결과</span></div><div><strong>전극 10개</strong><span>복합 fixture · 조합 16개</span></div><div><strong>13.48 s</strong><span>1000² 복합 short 기준 조건</span></div></div>'
    + table(["구분", "개발·검증 내용"], [
        ["구현", "입력 검증 → 2.5D 접촉/open → 전극 조합별 CA → 중복 제거 union CA → 결과 저장을 연결했다."],
        ["핵심 개선", "기존 2.5D 계산에 없던 xy와 높이 차의 <b>대각 거리 CA</b>를 추가하고 공정 산포 계산 경로를 마련했다."],
        ["복합 검증", "한 layer의 여러 전극, 다층 정렬·부분 중첩·대각 교차, 음성 대조와 hotspot 중첩을 확인했다."],
        ["결과 해석", "현재 값은 기하 접촉 위험이다. 실제 불량률·수율 또는 제품 승인값으로 해석하지 않는다."],
        ["다음 순서", "기존 sim 실제 동일 입력 비교 → 대체 판정 → 실제 제품 산포/Fab → particle-sim-compare 연동."],
    ])
    + source("README·SPEC·현재 구현 / evidence/pytest-p0.log, complex-validation.json, benchmark-*.json."),
)

page(
    "02 / 개선 목표",
    "기존 2.5D 계산에 대각 거리와 산포를 보완한다",
    '<p class="lead">기존 방식도 전극 높이를 반영한 2.5D 계산이었다. 부족했던 부분은 이물 중심에서 전극까지의 xy 거리와 높이 차를 하나의 대각 거리로 결합하지 못했고, mask 공정 산포를 반영하지 못한 점이다.</p>'
    + table(["기존 방식", "제안 방식", "의미"], [
        ["높이 반영, xy·z 분리 판단", "xy와 z를 결합한 대각 거리", "기존 계산 대상이 아니던 대각 CA를 추가"],
        ["고정 형상 한 조건", "mask별 overlay·CD 조건 반복", "위치별 조건부 hotspot과 평균 CA"],
        ["전극 조합별 결과", "조합별 결과와 OR union을 함께 보존", "원인 전극과 구조 전체 위험을 모두 확인"],
        ["Fab 분포와 계산 혼합 가능", "크기별 시뮬레이션 후 Fab 자료 적용", "구조 비교와 공정 분포 효과를 분리"],
    ])
    + '<div class="flow"><div>입력<br>mask·stack·전극</div><i>→</i><div>2.5D 계산<br>대각 거리</div><i>→</i><div>산포 반복<br>overlay·CD</div><i>→</i><div>사후 적용<br>Fab 이물 분포</div></div>'
    + '<h2>우선 비교할 값</h2><p>첫 판단은 동일한 이물 분포를 적용한 <b>구조 상대 비교</b>다. 실제 Fab 자료는 layer별·이물 크기별 평균 이물 개수와 검사면적을 사용하며, 크기별 CA 계산이 끝난 뒤 연결한다.</p>'
    + note("기존 결과가 없는 G–A 같은 대각 전극 조합은 기존값 0이 아니라 <b>N/A · new pair</b>로 표시한다.")
    + source("REPORT-SCRIPT 2페이지 확정 내용 / CODE-CHANGE-BACKLOG의 기존 sim 비교 계약."),
)

page(
    "03 / 입력과 출력",
    "적층 높이와 평가 전극 조합을 함께 읽는다",
    '<p class="lead">한 layer에는 여러 전극이 있을 수 있고 다른 layer의 전극과 복잡하게 겹칠 수 있다. 각 전극은 독립 mask로 입력하며, 전극 조합별 결과와 모든 조합을 중복 제거한 union 결과를 함께 저장한다.</p>'
    + fig("explainer_layout_3d.png", "그림 1. 왼쪽은 mask1·mask2·mask3의 상하 배치, 오른쪽은 평가 전극 조합이다. A–B는 같은 mask2의 lateral short, G–A와 B–P는 서로 다른 높이의 vertical short다.", "wide layout3d")
    + table(["입력", "쉽게 읽는 의미", "주요 출력"], [
        ["제품 manifest + 이진 mask", "좌표·축척, 적층, 전극 ID·net, 이물 유입 공정", "검증된 계산 형상"],
        ["이물 크기 + 평가 전극 조합", "어떤 크기의 이물이 어느 전극 사이를 연결하는지", "조합별 hotspot·CA"],
        ["mask별 산포 YAML", "실측 overlay x/y와 CD 분포", "조건부 확률 지도·평균 CA"],
        ["Fab 이물 표", "layer·크기별 평균 개수와 검사면적", "layer·전극 조합별 가중 위험"],
    ])
    + note("전극 mask의 합집합은 layer mask의 <b>부분집합</b>일 수 있다. 남는 TR 영역은 기록하되, 신호 분리 기능을 이물 불량 판정에 넣을지는 후속 물리 검토 항목으로 둔다.")
    + source("evidence/case/manifest.yaml, masks/ / figures/explainer_layout_3d.png."),
)

page(
    "04 / 계산 방법 · 접촉",
    "두 전극까지의 2.5D 대각 거리를 각각 계산한다",
    '<p><b>판정 순서:</b> ① 이물 직경에서 반지름과 중심 높이를 정한다. ② 중심에서 각 전극까지의 평면 거리와 높이 차를 구한다. ③ 두 거리를 결합한 대각 거리 Q를 계산한다. ④ 두 전극 모두 접촉하면 해당 중심 위치를 hotspot으로 기록한다.</p>'
    + '<div class="formula">r = d / 2 &nbsp; · &nbsp; Q = √(D<sub>xy</sub>² + D<sub>z</sub>²)<br>margin = r + contact margin − Q &nbsp; · &nbsp; hotspot = (margin<sub>A</sub> ≥ 0) AND (margin<sub>B</sub> ≥ 0)</div>'
    + fig("explainer_contact.png", "그림 2. mask1 형성 후 유입된 d=10 µm 이물의 G–A 단면 예시. x축은 이물 중심 위치, y축은 두 전극 중 더 작은 접촉 여유다. 0 이상인 주황 구간만 두 전극 동시 접촉 위치다.", "wide")
    + table(["구분", "해석"], [
        ["lateral short", "같은 layer의 서로 다른 net 전극 두 개에 동시에 닿는 위치"],
        ["vertical short", "서로 다른 layer의 전극 두 개에 동시에 닿는 위치. 상하 높이 확인이 필수"],
        ["contact margin", "첫 비교에서는 0 µm. 재질별 전기 접촉 거리 차이는 후속 검토"],
        ["이물 형상", "현재 원형/구형만 사용. 비원형 이물은 후속 검토"],
    ])
    + note("vertical short는 기하학상 접촉 가능한 위치다. 전도성, 절연막 파괴와 접촉저항은 아직 계산하지 않는다.")
    + source("geometry/contact.py, particles.py, surfaces.py / failure/lateral_short.py, vertical_short.py."),
)

page(
    "05 / 계산 방법 · open",
    "이물이 배선을 가려 연결을 끊는 위치를 찾는다",
    '<p class="lead">open은 배선끼리 접촉해서 생기는 불량이 아니다. 이물이 배선을 가리는 크기와 위치 때문에, 남은 도체에서 반드시 연결되어야 할 두 단자 영역 사이의 경로가 사라지는 불량이다.</p>'
    + '<div class="flow"><div>배선 + 두 단자 영역</div><i>→</i><div>후보 중심에 이물 배치</div><i>→</i><div>가려진 도체 제거</div><i>→</i><div>두 영역의 연결 검사</div></div>'
    + fig("explainer_open.png", "그림 3. 직선·우회·분기·ㄱ자 배선의 open 전수 대조. 주황색은 손상된 배선이 아니라, 그 위치에 이물 중심이 오면 open이 되는 후보 위치다.", "open")
    + '<p><b>0 px의 의미:</b> hotspot이 0이라는 뜻이 아니다. 구현 결과와 24×24 격자의 모든 중심 위치를 직접 계산한 기준 결과 사이의 불일치가 0픽셀이라는 뜻이다.</p>'
    + table(["입력 방식", "현재/계획", "이유"], [
        ["두 점 좌표", "합성 시험·이전 입력 호환", "실제 단자 형상을 한 픽셀로 대표하기 어려움"],
        ["두 단자 영역 mask", "실제품 기본 입력 · 구현 완료", "영역 일부가 남아도 경로가 존재하는지 판정 가능"],
        ["ROI 경계 fallback", "open proxy로만 유지", "실제 전기 단자를 지정한 결과와 구분"],
    ])
    + source("failure/open_failure.py / tests/unit/test_open_failure.py / evidence/recheck.json."),
)

page(
    "06 / 계산 방법 · 산포",
    "산포 조건마다 hotspot을 다시 계산해 위치별 비율을 구한다",
    '<p>첫 실제 적용에서는 Fab 실측값이 있는 <b>mask별 overlay x·y와 CD</b>만 산포로 사용한다. CD와 overlay는 서로 독립이며 overlay x·y도 독립으로 처리한다. 근거가 없는 두께·이물 크기·중심 이동 산포는 0으로 고정한다.</p>'
    + '<div class="formula">P<sub>geom</sub>(x,y | layer,d,pair) = hotspot이 된 산포 조건 수 / 전체 조건 수<br>CA<sub>mean</sub> = Σ P<sub>geom</sub>(x,y) × 픽셀 면적</div>'
    + fig("explainer_probability_profile.png", "그림 4. mask1 형성 후 d=10 µm 조건의 위치별 접촉 비율. 색 선은 전극 조합별 결과, 굵은 남색 선과 면은 각 산포 조건에서 조합들을 먼저 OR한 union 결과다.", "wide")
    + '<h2>Sobol 표본의 의미</h2><p>Sobol은 0~1 누적확률 공간을 고르게 채운 뒤 각 값을 지정 분포의 역누적분포로 변환한다. 정규분포를 지정하면 최종 overlay·CD 변화량은 정규분포를 근사한다. Sobol 자체가 균일하다는 이유로 실제 변화량이 균일분포가 되는 것은 아니다.</p>'
    + note("수렴은 중첩 표본 64→128→256으로 확인하고 필요하면 512로 늘린다. 최종 조건 수에서는 scramble seed를 3개 이상 바꾸어 결과 변동도 확인한다.")
    + source("uncertainty / aggregation/probability / evidence/prob / CODE-CHANGE-BACKLOG 산포 입력 항목."),
)

page(
    "07 / Fab 자료 적용",
    "시뮬레이션 후 layer·크기별 이물 분포를 연결한다",
    '<p class="lead">Fab 자료는 layer별·이물 크기별 <b>평균 이물 개수</b>와 sample당 검사면적이다. 평균 개수를 검사면적으로 나누어 평균 밀도를 만든 뒤, 같은 layer와 크기의 CA에 적용한다.</p>'
    + '<div class="formula">평균 이물 밀도 [개/mm²] = 평균 이물 개수 / sample당 검사면적 [mm²]<br>가중 기하 위험 = CA [mm²] × 평균 이물 밀도 [개/mm²]</div>'
    + table(["관점", "계산", "답하는 질문"], [
        ["layer별", "해당 유입 공정의 모든 hotspot을 OR한 union CA × 밀도", "어느 Fab layer의 전체 위험이 큰가?"],
        ["전극 조합별", "pair CA × 밀도, 중복 배분 기여도 병행", "어느 전극 조합이 민감하고 주요한가?"],
        ["layer × 전극", "Fab layer × 이물 크기 × 실제 pair ID", "어느 공정 이물이 어느 전극 불량에 크게 연결되는가?"],
    ])
    + '<div class="cols">'
    + fig("probability_fab_comparison.png", "그림 5. 합성 Fab 분포를 적용한 기능 경로 예시. 현재 값은 pair 단순합이므로 최종 전체 위험값으로 사용하지 않는다.", "half")
    + '<div><h2>중복 처리</h2><p><b>raw pair risk</b>는 조합 자체의 민감도다. 한 위치에서 k개 조합이 동시에 hotspot이면 <b>배분 기여도</b>는 각 조합에 1/k씩 나눈다. 배분 기여도 합은 layer union과 일치한다.</p><p>이 배분은 물리적 원인을 확정하는 식이 아니라 중복 없이 설명하기 위한 보고 규칙이다.</p></div></div>'
    + note("가중값은 기하 기반 기대 위험이다. 실제 전기 불량률·예상 불량 개수·수율로 확정하려면 제품 대조가 필요하다.")
    + source("io/fab.py, aggregation/fab_weighting.py / second-test-and-io-guide §9.5."),
)

page(
    "08 / 합성 검증",
    "통과한 계약과 아직 확인하지 않은 범위를 구분한다",
    '<div class="metrics"><div><strong>279 passed</strong><span>전체 pytest 실행</span></div><div><strong>0 px</strong><span>open 전수 기준과 불일치</span></div><div><strong>0 µm²</strong><span>복합 CA 적분 최대 오차</span></div></div>'
    + table(["검증 항목", "확인 결과", "남은 범위"], [
        ["산포 0 → 기준 조건", "불량 단위·union hotspot과 CA 일치", "실제 제품 입력"],
        ["같은 seed·설정", "확률·union 배열 재현, 여러 seed 기록", "실제 Fab 값 수렴"],
        ["open", "단자 영역·적용 공정과 전수 대조 0 px", "실제 단자 mask·대형 ROI"],
        ["크기 단조성", "2→4→6→8→10 µm 위반 없음", "비원형 이물"],
        ["대각 위치 이동", "상부 전극 이동 방향으로 hotspot 이동", "기존 sim 실측 대조"],
        ["중복 제거", "OR·동시 단위 수와 1/k 배분 합계 일치", "실제 제품 위험 순위"],
    ])
    + '<h2>검증 단계의 판단</h2><p>합성 결과는 계산 규칙과 다중 전극 처리의 기능 근거다. 실제 제품의 유효성이나 기존 현장 sim 대체 가능성을 증명하지는 않는다. 다음 검증부터 공통 pair와 신규 대각 pair를 구분해 비교한다.</p>'
    + note("외부 기존 sim 결과와 내부 <code>legacy_failure</code>는 다르다. 내부 값은 xy kernel 참고 결과이며 기존 sim 대체 판정에 직접 사용하지 않는다.")
    + source("evidence/pytest-p0.log, recheck.json, complex-validation.json, regression-probes.json."),
)

triple_s2_d6 = triple[(triple.source_step == "S_M2") & (triple.particle_size_um == 6) & (triple.failure_mode == "any")].iloc[0]
page(
    "09 / CA와 중복 제거",
    "전극 조합별 원인과 구조 전체 위험을 함께 본다",
    '<h2>pilot 기준 조건 · d=10 µm</h2>'
    + table(["이물 유입 시점", "A–B lateral", "G–A vertical", "B–P vertical", "전체 union"], [
        ["mask1 형성 후", "700", "350", "700", "1,750 µm²"],
        ["mask2 형성 후", "100", "0", "400", "500 µm²"],
        ["평탄화막 형성 후", "0", "0", "0", "0 µm²"],
    ])
    + '<div class="cols">'
    + fig("nominal_ca_curve_S_M1.png", "그림 6. mask1 형성 후 이물 크기별 pair CA. 크기가 커질수록 접촉 가능한 중심 위치가 넓어진다.", "half")
    + fig("triple_union_vs_sum_S_M2.png", "그림 7. 세 전극 조합의 단순합과 union 비교. 같은 중심 위치의 중복은 union에서 한 번만 센다.", "half")
    + "</div>"
    + table(["겹침 예시", "union CA", "pair CA 단순합", "해석"], [
        ["mask2 형성 후 / 10 µm", "350 µm²", "650 µm²", "여러 조합이 같은 위치를 공유"],
        ["mask2 형성 후 / 6 µm", f"{float(triple_s2_d6.ca_union_nominal_um2):g} µm²", f"{float(triple_s2_d6.ca_pair_sum_nominal_um2):g} µm²", "겹침이 없으면 두 값이 같음"],
    ])
    + note("union CA는 해당 layer의 전극 mask 합집합이 아니다. 같은 이물 유입 공정·크기에서 선언된 lateral short, vertical short, open hotspot의 OR다.")
    + source("evidence/run/ca.parquet, ca_union.parquet / triple/ca_union.parquet."),
)

page(
    "10 / 2차 합성 시험",
    "다중 전극과 다층 중첩을 한 입력에서 검증했다",
    '<div class="cols">'
    + fig("complex_layout_overlay.png", "그림 8. mask1 3개, mask2 4개, mask3 3개 전극의 굴곡·분기·섬·교차 형상.", "half")
    + fig("complex_pair_S_M1_G2__D3_d10.png", "그림 9. mask1 형성 후 G2–D3, d=10 µm의 기준 2.5D hotspot과 접촉 여유.", "half")
    + "</div>"
    + table(["입력 항목", "2차 시험"], [
        ["격자", "200×200 px, 0.5 µm/px, ROI 100×100 µm"],
        ["전극", "10개 — mask1 3, mask2 4, mask3 3"],
        ["전극 조합", "16개 — lateral short 6, vertical short 10"],
        ["형상", "L·T 분기, 좁은 간격, 분리된 섬, 정렬·부분 중첩·대각 교차"],
        ["조건", "이물 유입 공정 3 × 직경 5개, open·산포 제외"],
    ])
    + '<h2>자동 대조 결과</h2>'
    + table(["검사", "결과"], [
        ["지정 vertical short 양성", "8개 조합 모두 hotspot 존재"],
        ["먼 거리 음성 대조", "G1–D4, G3–D3 모두 0 px"],
        ["크기 단조성 / union 일치", "위반 0 / OR 불일치 0"],
        ["최대 동시 영향 조합", str(complex_validation["max_simultaneous_pairs"])],
        ["내부 kernel 차이", "2,911 px, 미분류 0"],
    ])
    + note("복합 fixture는 기하학상 가능한 예시다. 실제 제품 mask와 기존 sim 결과를 대조해야 대체 판단에 사용할 수 있다.")
    + source("evidence/complex-case, complex-run, complex-validation.json."),
)

page(
    "11 / 2차 시험 · TAT",
    "형상 복잡도와 이미지 크기의 비용을 분리했다",
    '<p>동일 장비에서 각 조건을 1회 워밍업 후 5회 실행했다. 시간은 wall time 중앙값, 메모리는 5회 중 최대 RSS다. 이물 유입 공정 3개·직경 5개·open 제외 조건으로 비교했다.</p>'
    + fig("complex_performance.png", "그림 10. 기준 조건 TAT와 최대 메모리. 색은 형상 종류, x축은 pilot/복합과 200²/1000²의 네 조건만 표시했다.", "wide performance")
    + table(["입력", "전극 / 조합", "TAT", "최대 RSS", "출력"], [
        ["pilot 200²", "4 / 3", f'{bench["pilot-200"]["wall_median_seconds"]:.2f} s', f'{bench["pilot-200"]["peak_rss_max_bytes"]/1e6:.1f} MB', f'{bench["pilot-200"]["output_bytes_median"]/1e6:.2f} MB'],
        ["복합 200²", "10 / 16", f'{bench["complex-200"]["wall_median_seconds"]:.2f} s', f'{bench["complex-200"]["peak_rss_max_bytes"]/1e6:.1f} MB', f'{bench["complex-200"]["output_bytes_median"]/1e6:.2f} MB'],
        ["pilot 1000²", "4 / 3", f'{bench["pilot-1000"]["wall_median_seconds"]:.2f} s', f'{bench["pilot-1000"]["peak_rss_max_bytes"]/1e6:.1f} MB', f'{bench["pilot-1000"]["output_bytes_median"]/1e6:.2f} MB'],
        ["복합 1000²", "10 / 16", f'{bench["complex-1000"]["wall_median_seconds"]:.2f} s', f'{bench["complex-1000"]["peak_rss_max_bytes"]/1e6:.1f} MB', f'{bench["complex-1000"]["output_bytes_median"]/1e6:.2f} MB'],
    ])
    + f'<p class="small">복합 1000²는 pilot 1000²보다 TAT 약 {bench["complex-1000"]["wall_median_seconds"]/bench["pilot-1000"]["wall_median_seconds"]:.2f}배, 최대 RSS 약 {bench["complex-1000"]["peak_rss_max_bytes"]/bench["pilot-1000"]["peak_rss_max_bytes"]:.2f}배다. 출력 크기는 전극 조합 수와 Zarr 압축률의 영향을 함께 받는다.</p>'
    + source("evidence/benchmark-*.json / benchmark_nominal.py / make_complex_figures.py."),
)

page(
    "12 / 대형·산포 성능",
    "기준 조건과 산포 계산의 TAT를 따로 읽는다",
    '<p class="lead">장비는 Mac mini M4 Pro, CPU 12코어, 메모리 48 GB다. 아래 값은 단일 대형 합성 입력의 실측이며 실제 제품 최대치가 아니다.</p>'
    + table(["실행 조건", "TAT", "최대 RSS", "해석"], [
        ["1000×1000, 직경 32개, 공정 3, pair 3<br>기준 short · open 제외", f'{large_short["wall_seconds"]:.2f} s', f'{large_short["peak_rss_bytes"]/1e6:.1f} MB', "크기 수 증가 기준선"],
        ["1000×1000, 직경 5개, 64 산포 조건<br>short · open 제외", f'{large_prob["wall_seconds"]:.2f} s', f'{large_prob["peak_rss_bytes"]/1e6:.1f} MB', "64조건 전체 반복"],
        ["같은 산포 입력 · P0 메모리 보완 후", f'{large_prob_p0["wall_seconds"]:.2f} s', f'{large_prob_p0["peak_rss_bytes"]/1e6:.1f} MB', "TAT 약 +5%, RSS 약 −41%"],
    ])
    + '<h2>200×200 복합 형상의 failure mode별 파일럿</h2>'
    + table(["실행", "활성 조건", "TAT 최초→최종", "최종 TAT 비율", "통합 CA 비율"], [
        ["lateral short", "90", f'{mode_tat["lateral_short"]["wall_median_seconds"]:.2f}→{mode_tat_final["lateral_short"]["wall_median_seconds"]:.2f} s', f'{mode_tat_final["lateral_short"]["isolated_tat_share"]*100:.2f}%', f'{mode_ca["lateral short"]["allocated_ca_share"]*100:.2f}%'],
        ["vertical short", "150", f'{mode_tat["vertical_short"]["wall_median_seconds"]:.2f}→{mode_tat_final["vertical_short"]["wall_median_seconds"]:.2f} s', f'{mode_tat_final["vertical_short"]["isolated_tat_share"]*100:.2f}%', f'{mode_ca["vertical short"]["allocated_ca_share"]*100:.2f}%'],
        ["open", "35", f'{mode_tat["open"]["wall_median_seconds"]:.2f}→{mode_tat_final["open"]["wall_median_seconds"]:.2f} s', f'{mode_tat_final["open"]["isolated_tat_share"]*100:.2f}%', f'{mode_ca["open"]["allocated_ca_share"]*100:.2f}%'],
        ["통합", "275", f'{mode_tat["combined"]["wall_median_seconds"]:.2f}→{mode_tat_final["combined"]["wall_median_seconds"]:.2f} s', "4.04배 향상", "100.00%"],
    ])
    + note("공정 간 open 결과 재사용과 전극 bounding box 연결성 판정으로 open은 5.27배, 통합은 4.04배 빨라졌다. 최적화 전후 hotspot·CA는 픽셀 단위로 같다. CA 비율은 Fab 가중 전 합성 기하 결과다.")
    + f'<p class="small"><b>동일 형상 padded 1000²:</b> lateral {padded_1000["lateral_short"]["wall_median_seconds"]:.2f}초 · vertical {padded_1000["vertical_short"]["wall_median_seconds"]:.2f}초 · open {padded_1000["open"]["wall_median_seconds"]:.2f}초 · 통합 {padded_1000["combined"]["wall_median_seconds"]:.2f}초 / 최대 RSS {padded_1000["combined"]["peak_rss_max_bytes"]/1e6:.1f}MB / 출력 {padded_1000["combined"]["output_bytes_median"]/1e6:.1f}MB. 같은 국부 형상을 큰 빈 ROI 중앙에 둔 확장성 시험이며 실제 긴 배선 대표값은 아니다. 출력 대부분은 pair_margin_um이다.</p>'
    + source("evidence/large-short.json, large-prob*.json, benchmark-failure-modes-200*.json, benchmark-failure-modes-padded-1000.json, open-optimization-equivalence.json / benchmark_failure_modes.py."),
)

page(
    "13 / 회귀 검토",
    "작은 반례와 전수 대조로 계산 경계를 확인했다",
    '<p>계약 검증만으로 기하 계산의 정확성을 보장할 수 없다. 작은 합성 반례에서 모든 픽셀을 직접 대조하고, 발견한 오류를 회귀 테스트로 고정했다.</p>'
    + table(["대표 결함", "발생 가능한 오판", "현재 조치"], [
        ["중심 offset x/y 뒤바뀜", "비등방 격자에서 다른 방향으로 hotspot 이동", "축 순서 수정·회귀 테스트"],
        ["ㄱ자 open 후보 누락", "실제로 단절되는 위치를 정상으로 누락", "모서리 후보 보완·전수 비교"],
        ["비등방 CD 반픽셀", "x/y 물리 길이가 다른 경계를 같은 길이로 보정", "축별 픽셀 길이 적용"],
        ["Fab bin 행 순서", "같은 크기가 다른 밀도 구간에 연결", "정렬·반개구간 규칙 고정"],
        ["union 수렴 누락", "pair 수렴만으로 전체 위험 수렴 오인", "union 계열 수렴 기록 추가"],
    ])
    + '<h2>P0 코드 보완 회귀 결과</h2>'
    + table(["보완", "확인 결과"], [
        ["layer mask 부분집합", "전극 밖 TR 영역 허용·면적 기록, layer 밖 전극은 오류"],
        ["open 단자 영역", "점 좌표 호환, 영역 mask 연결성과 적용 공정 검증"],
        ["배분 기여도", "비중첩·완전 중첩·부분 중첩에서 합계가 union과 일치"],
        ["산포 union 열", "기준 CA와 산포 평균·분위수를 다른 열에 보존"],
    ])
    + source("tests/unit / evidence/recheck.json, regression-probes.json / CODE-CHANGE-BACKLOG.md."),
)

page(
    "14 / 폐쇄망 재현",
    "입력·코드·출력의 버전을 함께 고정한다",
    '<p class="lead">현재 저장소에는 <code>pyproject.toml</code>과 <code>uv.lock</code>이 있다. 폐쇄망 반입은 대상 OS·CPU·Python용 패키지 묶음과 빈 환경 재설치 시험이 추가로 필요하다.</p>'
    + table(["반입 항목", "현재 확인", "필요 조치"], [
        ["Python", "3.12.13 실행 근거", "대상 환경 버전 고정"],
        ["의존성", "uv.lock 존재", "대상 플랫폼 wheelhouse 준비"],
        ["신규 설치", "미측정", "빈 환경 설치 시간·로그·hash 보존"],
        ["회귀 테스트", "279개 통과", "실제 비교 adapter 추가 시 재실행"],
        ["실행 추적", "입력·출력 hash, git commit/dirty, build ID 기록", "clean commit 실행 보존"],
    ])
    + '''<h2>반입 후 기본 확인 절차</h2><pre>python3.12 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheelhouse -r requirements.lock.txt
export PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1
.venv/bin/python -m pytest -q
.venv/bin/python -m particle_sim validate-input --product-config product.yaml</pre>'''
    + '<p>최종 비교 manifest에는 기존 sim과 particle-sim의 버전, 입력·출력 hash와 좌표 변환 규칙을 함께 기록한다.</p>'
    + source("pyproject.toml, uv.lock / evidence/environment.json, run/run_manifest.json."),
)

page(
    "15 / 우선순위",
    "P0 보완을 마치고 기존 sim 실제 대조로 이동한다",
    table(["순서", "작업", "완료 기준"], [
        ["1 · 완료", "입력·출력 계약 보완", "TR 잔여 영역, 조합 감사, open 단자, 산포/union, Fab 연결"],
        ["2 · 완료", "회귀 검증·비교 엔진", "279개 테스트, schema·commit·출력 hash, compare-sim round-trip"],
        ["3 · 다음", "기존 sim 실제 대조", "공통 조합 hotspot·XOR·CA 차이 분류, 신규 대각 조합은 N/A"],
        ["4 · P1", "실제품 기준·산포 pilot", "mask별 실측 overlay·CD, 64→128→256 수렴과 여러 seed"],
        ["5 · P1", "Fab 자료 사후 적용", "평균 개수/검사면적을 layer·크기별 CA에 연결"],
        ["6 · P1", "particle-sim-compare", "layer·전극·교차 위험 순위와 수기 계산 일치"],
    ])
    + '<h2>추후 물리 검토 항목</h2>'
    + table(["항목", "현재 처리", "추후 판단"], [
        ["TR 신호 분리", "기하 mask 잔여 영역만 기록", "불량 판정에 전기 기능을 넣을지 결정"],
        ["전도성/비전도성", "contact margin 0", "재질별 접촉 거리·불량 모드"],
        ["이물 형상", "원형/구형", "타원·섬유·불규칙 형상"],
        ["예상 불량", "기하 위험까지만", "제품 불량 자료로 변환식 검증"],
    ])
    + note("P0 구현 상태와 실제 비교에 남은 조건은 <code>docs/CODE-IMPLEMENTATION-PLAN.md</code>에 기록했다.")
    + source("REVIEW-PLAN.md / CODE-CHANGE-BACKLOG.md."),
)

manifest_left = """coordinate_system:
  pixel_size_x_um: 0.5
  pixel_size_y_um: 0.5
  roi_width_um: 100.0
  roi_height_um: 100.0

particle_sizes_um: [2, 4, 6, 8, 10]
embedding_depth_um: 0.0
contact_margin_um: 0.0"""
manifest_right = """particle_source_steps:
- step_id: after_mask1
  lands_after_layer: mask1
  fab_distribution_layer: fab_mask1

electrode_pairs:
- pair_id: gate1__data1
  electrode_a: gate1
  electrode_b: data1
  failure_mode: vertical_short"""
page(
    "16 / 입력 가이드",
    "제품 형상·산포·Fab 자료를 서로 다른 입력으로 관리한다",
    '<p>아래는 필드 의미를 보여주는 발췌이며 실행 가능한 전체 YAML은 아니다. 입력 shape는 고정되지 않으며 모든 mask의 크기와 좌표계가 같으면 1000×1000 이상도 사용할 수 있다.</p>'
    + f'<div class="cols"><pre>{html.escape(manifest_left)}</pre><pre>{html.escape(manifest_right)}</pre></div>'
    + table(["입력 묶음", "담는 내용", "검증 관계"], [
        ["제품 manifest + mask", "좌표·적층·전극·net·평가 조합", "전극은 layer mask의 부분집합, 같은 layer 전극끼리 비중첩"],
        ["산포 YAML", "mask별 overlay x/y·CD 분포와 근거", "첫 실제 계산은 실측값만, 변수들은 독립"],
        ["Fab 이물 표", "layer·크기 bin별 평균 개수·검사면적·기간", "시뮬레이션 후 평균 밀도로 변환"],
        ["open 단자", "시작·도착 단자 영역 mask", "비어 있지 않고 배선 내부이며 서로 분리"],
    ])
    + '<h2>입력에서 반드시 구분할 이름</h2><p><b>mask1·mask2</b>는 전극이 놓인 물리 layer이고, <b>after_mask1</b>은 이물이 들어온 시점이며, <b>fab_mask1</b>은 Fab 이물 자료가 측정된 공정 이름이다. 세 값을 같은 문자열이라고 가정하지 않는다.</p>'
    + note("YAML은 # 주석을 허용하지만 주석은 실행 hash와 결과에 전달되지 않는다. 필드 의미와 근거는 이 가이드와 run manifest에서 관리한다.")
    + source("second-test-and-io-guide.md / evidence/case/manifest.yaml / schema.py."),
)

page(
    "17 / 출력 가이드",
    "layer·전극·교차 위험을 모두 확인할 수 있어야 한다",
    '<div class="cols"><div><h2>Zarr · 위치 지도 묶음</h2><p>조건별 2차원 hotspot 지도를 여러 장 쌓아 둔 형식이다. 큰 지도를 chunk로 나누므로 필요한 공정·크기·전극 조합과 위치만 읽는다.</p><p><b>질문:</b> 어디가 위험하며 여러 조합이 어디서 겹치는가?</p></div><div><h2>Parquet · 분석용 표</h2><p>행과 열이 있는 압축 표다. Fab layer, 이물 크기, 전극 ID, CA와 단위를 정렬·필터·집계한다. 사람이 볼 상위 순위는 CSV로도 만든다.</p><p><b>질문:</b> 어느 조건·layer·전극의 값이 큰가?</p></div></div>'
    + table(["최종 출력", "핵심 내용", "용도"], [
        ["hotspots.zarr / probability.zarr", "pair·union의 위치별 0/1 또는 조건부 비율", "공간 근거 확인"],
        ["ca.parquet / ca_union.parquet", "pair CA와 중복 제거 union CA", "기하 결과·수렴 확인"],
        ["risk_by_layer.parquet", "Fab layer별 union CA × 이물 밀도", "layer 위험 순위"],
        ["risk_by_pair.parquet", "Fab layer × 실제 전극 pair, raw·배분·고유 기여", "전극 조합·교차 위험 순위"],
        ["risk_by_open_path.parquet", "open 경로별 raw·배분·고유 기여", "주요 open 경로"],
        ["risk_by_electrode.parquet", "개별 전극 관여 점수와 합계 일치 기여", "반복 등장하는 주요 전극"],
        ["risk_summary.csv", "layer·조합·open·전극 상위 항목", "사람이 바로 여는 검토표"],
        ["run/risk manifest", "입력·코드·schema·출력 hash와 가정", "실행 재현"],
    ])
    + '<h2>현재 출력과 연동 보완</h2><p>현재 <code>ca.parquet</code>에는 실제 전극 ID가 있어 pair 분석의 기반은 있다. 다만 particle-sim-compare 내보내기는 역할·layer 중심이고, 실측 평균 이물 개수 가중과 산포 probability 내보내기가 아직 완성되지 않았다. 연동 시 실제 전극 ID·pair ID·Fab layer·배분 기여도를 구조화된 필드로 추가한다.</p>'
    + note("raw pair risk는 조합 민감도이므로 합산하지 않는다. 배분 기여도는 같은 위치의 k개 조합에 1/k를 적용해 합계가 layer union과 일치하도록 한다.")
    + source("second-test-and-io-guide §9 / CODE-CHANGE-BACKLOG §6–7 / particle-sim-compare data contract."),
)

base_css = """@page { size: A4; margin: 18mm 16mm; }
*{box-sizing:border-box}body{margin:0;color:#193247;background:#fff;font-family:"Apple SD Gothic Neo",AppleGothic,sans-serif;font-size:10pt;line-height:1.5;-webkit-print-color-adjust:exact;print-color-adjust:exact}.page{height:260mm;position:relative;page-break-after:always;padding:0 0 12mm}.page:last-child{page-break-after:auto}header{font-size:8pt;letter-spacing:.8px;color:#397c88;border-bottom:1px solid #b7cdd3;padding-bottom:3mm;margin-bottom:6mm;font-weight:700}header span{float:right;color:#6c7f8b;font-weight:400;letter-spacing:0}.kicker{font-size:9pt;color:#007f88;font-weight:700;letter-spacing:1px}h1{font-size:22pt;letter-spacing:-.7px;line-height:1.25;margin:2mm 0 5mm;color:#133148}h2{font-size:11.5pt;line-height:1.3;margin:4mm 0 2mm;color:#007c83}p{margin:2.2mm 0 3mm;word-break:keep-all}.lead{font-size:11.5pt;line-height:1.58}.small{font-size:8.5pt}table{width:100%;border-collapse:collapse;table-layout:fixed;margin:2.5mm 0 3mm;font-size:8.8pt;page-break-inside:avoid;break-inside:avoid}th{text-align:left;background:#193d53;color:white;font-weight:600}td,th{padding:1.9mm 2.1mm;border-bottom:1px solid #dce5e9;vertical-align:top;word-break:keep-all;overflow-wrap:anywhere}tr:nth-child(even) td{background:#f1f6f7}table th:first-child,table td:first-child{width:25%}.note{background:#eaf5f3;border-left:3px solid #159287;padding:2.7mm 3.5mm;margin:3mm 0;font-size:9.4pt;break-inside:avoid}.metrics{display:flex;gap:4mm;margin:5mm 0}.metrics>div{flex:1;border-top:3px solid #10a394;background:#f1f7f8;padding:3.5mm}.metrics strong{display:block;font-size:21pt;line-height:1.2;color:#153d53}.metrics span{font-size:8.8pt}.flow{display:flex;gap:2mm;align-items:center;margin:4mm 0}.flow div{flex:1;background:#edf4f6;padding:2.8mm;text-align:center;border:1px solid #c4d8de;font-size:8.8pt}.flow i{font-style:normal;color:#009b91}.cols{display:flex;gap:5mm;align-items:flex-start}.cols>*{width:calc(50% - 2.5mm);min-width:0}.cols table{font-size:8.4pt}.cols td,.cols th{padding:1.7mm}figure{margin:2.5mm 0 3mm;break-inside:avoid;text-align:center}figure img{width:100%;object-fit:contain;display:block;max-height:88mm}figure.wide img{max-height:58mm}figure.half img{height:57mm}figure.open img{max-height:82mm}figure.layout3d img{max-height:69mm}figure.performance img{max-height:67mm}figcaption{font-size:8pt;color:#4e6677;line-height:1.35;margin-top:1.5mm;text-align:left}figcaption b{color:#527583;font-weight:500}.formula{background:#f0f5f8;border:1px solid #ccdce4;padding:3.5mm;text-align:center;font-size:11pt;line-height:1.55;margin:3mm 0}code{font-family:Menlo,monospace;font-size:.84em;overflow-wrap:anywhere}pre{white-space:pre-wrap;word-break:break-all;font:7.8pt/1.4 Menlo,monospace;background:#f0f4f6;padding:3mm;border:1px solid #d7e3e7;margin:2.5mm 0}.source{font-size:7.3pt;color:#6a7c89;margin-top:3mm;overflow-wrap:anywhere}footer{position:absolute;bottom:0;left:0;right:0;border-top:1px solid #cbdce2;padding-top:2.5mm;display:flex;justify-content:space-between;font-size:8pt;color:#71838c}footer b{color:#007d82}@media screen{body{background:#e9eff2}.page{background:white;width:210mm;height:297mm;margin:8mm auto;padding:18mm 16mm 30mm;box-shadow:0 3px 20px #bccbd3}footer{bottom:18mm;left:16mm;right:16mm}}"""
css = base_css + "\n" + (R / "report-layout.css").read_text(encoding="utf-8")
document = (
    '<!doctype html><html lang="ko"><head><meta charset="utf-8">'
    '<meta name="viewport" content="width=device-width, initial-scale=1">'
    '<title>particle-sim 개발 결과 보고서</title><style>'
    + css + "</style></head><body>" + "".join(pages) + "</body></html>"
)
(R / "report.html").write_text(document, encoding="utf-8")
print(f"pages {len(pages)}")
