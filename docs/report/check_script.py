"""발표 원고(REPORT-SCRIPT.md)의 구조와 수치를 보고서와 대조한다.

저장소 루트에서 실행한다.

    python3 docs/report/check_script.py

표준 라이브러리만 쓰므로 폐쇄망에서도 그대로 동작한다. 구조 검사에서 하나라도 실패하면
종료 코드 1이다. 수치 대조는 사람이 판단할 항목이므로 경고로만 보고한다.

원고를 고친 뒤에는 이 스크립트를 먼저 돌리고, 새로 뜨는 수치가 있으면 근거를 확인해
`VERIFIED_NUMBERS`에 이유와 함께 추가한다.
"""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path

SCRIPT = Path("docs/report/REPORT-SCRIPT.md")
REPORT = Path("docs/report/report.html")
SRC = Path("src/particle_sim")

PAGE_COUNT = 25

# 쪽마다 반드시 있어야 하는 절
REQUIRED_SECTIONS = ("핵심 메시지", "발표 스크립트", "사실과 한계", "첨삭 요청")

# 원고가 인용하는 코드 식별자. 소스에서 사라지면 원고가 거짓이 된다.
CODE_IDENTIFIERS = (
    "CNO_COVERAGE_RULE",
    "cno_coverage_rule",
    "risk_by_contact",
    "insufficient_terminals",
    "baseline_open",
    "terminal_lost",
    "contact_lost",
    "auto_open",
    "new_unit",
)

# report.html에 없지만 근거를 확인한 수치. 이유 없이 늘리지 않는다.
VERIFIED_NUMBERS = {
    "6.36": "102.25÷16.07. '이렇게 계산하면 안 된다'는 반례로만 등장한다",
    "1.0": "64조건 중 64회 접촉의 계산 예시",
    "0.0625": "64조건 중 4회 접촉(4/64)의 계산 예시",
    "0.25": "대형 실행의 픽셀 크기. evidence/large-case/manifest.yaml의 pixel_size_x_um",
    "250": "대형 실행의 ROI. evidence/large-case/manifest.yaml의 roi_width_um",
    "1.49": "padded 1000² open 4.67초 ÷ 200² open 3.14초. HANDOFF-2026-09-19.md 8절에 명시",
}

# 원고 자체를 설명하는 수·비율이라 보고서와 대조할 대상이 아니다.
META_SECTIONS = (
    "## 이 문서의 상태",
    "## 17쪽 판과의 대응",
    "## 발표 시간 배분 제안",
    "## 용어 표기 원칙",
    "## 발표 전 공통 주의",
    "## 보고서 재생성 절차",
)


def fail(results: list[tuple[str, bool, str]], name: str, detail: str) -> None:
    results.append((name, False, detail))


def ok(results: list[tuple[str, bool, str]], name: str, detail: str = "") -> None:
    results.append((name, True, detail))


def body_without_meta(text: str) -> str:
    """원고 설명용 머리말·꼬리말을 뺀 본문만 남긴다."""
    out, skip = [], False
    for line in text.split("\n"):
        if line.startswith("## "):
            skip = any(line.startswith(m) for m in META_SECTIONS)
        if not skip:
            out.append(line)
    return "\n".join(out)


def check_structure(text: str, results: list) -> None:
    pages = re.findall(r"^## (\d+)페이지", text, re.M)
    expected = [str(i) for i in range(1, PAGE_COUNT + 1)]
    if pages == expected:
        ok(results, "쪽 구성", f"{PAGE_COUNT}개, 1~{PAGE_COUNT} 연속")
    else:
        fail(results, "쪽 구성", f"{len(pages)}개, 순서 {pages[:5]}…")

    missing = []
    for n in expected:
        m = re.search(rf"^## {n}페이지[^\n]*$(.*?)(?=^## |\Z)", text, re.M | re.S)
        if not m:
            missing.append(f"{n}쪽 없음")
            continue
        heads = re.findall(r"^### ([^\n]+)$", m.group(1), re.M)
        for want in REQUIRED_SECTIONS:
            if want not in heads:
                missing.append(f"{n}쪽 '{want}'")
    if missing:
        fail(results, "쪽별 필수 절", ", ".join(missing))
    else:
        ok(results, "쪽별 필수 절", f"{len(REQUIRED_SECTIONS)}개 × {PAGE_COUNT}쪽")


def check_appendix(text: str, results: list) -> None:
    m = re.search(r"^## 부록 [^\n]*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    if not m:
        ok(results, "부록", "없음 (선택 사항)")
        return
    empty = [
        h.group(1)
        for h in re.finditer(r"^### ([^\n]+)\n(.*?)(?=^### |\Z)", m.group(1), re.M | re.S)
        if not h.group(2).strip()
    ]
    if empty:
        fail(results, "부록", f"본문 없는 절: {empty}")
    else:
        n = len(re.findall(r"^### ", m.group(1), re.M))
        ok(results, "부록", f"{n}개 절, 빈 절 없음")


def check_tables(text: str, results: list) -> None:
    lines = text.split("\n")
    bad, count, i = [], 0, 0
    while i < len(lines):
        if lines[i].startswith("|"):
            j = i
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            count += 1
            widths = {line.count("|") for line in lines[i:j]}
            if len(widths) != 1:
                bad.append(i + 1)
            i = j
        else:
            i += 1
    if bad:
        fail(results, "표 열 수", f"불일치 블록 시작 줄: {bad}")
    else:
        ok(results, "표 열 수", f"{count}개 블록 일관")


def check_delimiters(text: str, results: list) -> None:
    if text.count("“") != text.count("”"):
        fail(results, "따옴표 짝", f"{text.count(chr(0x201c))} vs {text.count(chr(0x201d))}")
    elif text.count("```") % 2:
        fail(results, "코드펜스 짝", "홀수 개")
    else:
        ok(results, "따옴표·코드펜스", f"따옴표 {text.count(chr(0x201c))}쌍, 펜스 짝 맞음")


def check_page_refs(text: str, results: list) -> None:
    """옛 17쪽 판 기준의 'N페이지' 상호참조가 본문에 남았는지 본다."""
    stale = [
        (i + 1, line.strip())
        for i, line in enumerate(body_without_meta(text).split("\n"))
        if re.search(r"\d+페이지", line) and not line.startswith("## ")
    ]
    if stale:
        fail(results, "쪽 번호 상호참조", f"{len(stale)}건: {stale[:3]}")
    else:
        ok(results, "쪽 번호 상호참조", "본문에 'N페이지' 참조 없음 (모두 '…쪽')")


def check_identifiers(results: list) -> None:
    if not SRC.is_dir():
        fail(results, "코드 식별자", f"{SRC} 없음")
        return
    blob = "\n".join(
        p.read_text(encoding="utf-8", errors="ignore") for p in SRC.rglob("*.py")
    )
    missing = [name for name in CODE_IDENTIFIERS if name not in blob]
    if missing:
        fail(results, "코드 식별자", f"소스에 없음: {missing}")
    else:
        ok(results, "코드 식별자", f"{len(CODE_IDENTIFIERS)}개 모두 존재")


def check_numbers(text: str, results: list) -> list[str]:
    if not REPORT.is_file():
        fail(results, "수치 대조", f"{REPORT} 없음")
        return []
    plain = html.unescape(re.sub(r"<[^>]+>", " ", REPORT.read_text(encoding="utf-8")))
    plain_nc = plain.replace(",", "")
    tokens = {
        t
        for line in body_without_meta(text).split("\n")
        if not line.startswith("## ")
        for t in re.findall(r"\d[\d,]*\.?\d*", line)
    }
    unknown = sorted(
        t
        for t in tokens
        if t not in plain and t.replace(",", "") not in plain_nc and t not in VERIFIED_NUMBERS
    )
    if unknown:
        fail(results, "수치 대조", f"근거 미확인 {len(unknown)}개 (아래 목록)")
    else:
        ok(results, "수치 대조", f"미일치 {len(VERIFIED_NUMBERS)}개 전부 확인된 예외")
    return unknown


def main() -> int:
    if not SCRIPT.is_file():
        print(f"{SCRIPT} 을 찾을 수 없다. 저장소 루트에서 실행한다.")
        return 2
    text = SCRIPT.read_text(encoding="utf-8")

    results: list[tuple[str, bool, str]] = []
    check_structure(text, results)
    check_appendix(text, results)
    check_tables(text, results)
    check_delimiters(text, results)
    check_page_refs(text, results)
    check_identifiers(results)
    unknown = check_numbers(text, results)

    width = max(len(name) for name, _, _ in results)
    print(f"{SCRIPT} · {len(text.splitlines())}줄\n")
    for name, passed, detail in results:
        print(f"  {'통과' if passed else '실패'}  {name:<{width}}  {detail}")

    if unknown:
        print("\n보고서에 없고 확인되지 않은 수치:")
        for t in unknown:
            print(f"  {t}")
        print("\n근거를 확인한 뒤 check_script.py의 VERIFIED_NUMBERS에 이유와 함께 추가한다.")

    failed = [name for name, passed, _ in results if not passed]
    print(f"\n{len(results) - len(failed)}/{len(results)} 통과")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
