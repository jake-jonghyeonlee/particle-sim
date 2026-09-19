# SPEC: particle-sim

`particle-sim`은 display layout의 2.5D 이물 기하 접촉 위험을 계산한다. 입력은 layer별 패턴,
분리된 전극 mask, 적층 두께, 이물 유입 공정과 직경이다. 출력은 `lateral_short`,
`vertical_short`, `open` hotspot과 critical area(CA), 중복 제거 union, Fab 이물 밀도를 적용한
layer·전극 조합·open·개별 전극 위험도다. 이 값은 기하 접촉 위험이며 전기적 불량확률이 아니다.

## 입력 계약

- Python 배열 축은 `(y, x)`, 물리 좌표는 `(x_um, y_um, z_um)`, 길이는 µm, 면적은 µm²다.
- 모든 PNG/TIFF mask는 같은 shape의 이진 배열이어야 한다. resize, threshold, 빈 값 보정은 하지 않는다.
- 한 layer에 여러 전극을 둘 수 있다. 같은 layer의 전극은 겹치면 안 되고 각 전극은 layer mask의
  부분집합이어야 한다. 전극 합집합에 포함되지 않은 TR·비평가 영역은 허용하고 면적을 기록한다.
- 서로 다른 net의 모든 전극 조합은 schema 2.0에서 `electrode_pairs` 또는 이유가 있는
  `excluded_electrode_pairs` 중 하나에 있어야 한다. 미분류 조합이 있으면 실행을 거부한다.
- 실제 제품 open은 두 단자 영역 mask와 `applicable_particle_source_steps`로 정의한다. 이물이 한
  단자를 모두 가리거나 두 단자의 4-연결 경로를 끊으면 open이다.
- overlay/CD schema 2.0은 x/y 분포, 기준 layer, 부호, 항목별 실측/가정 근거를 기록한다. 현재 변수는
  독립으로 표본화한다.
- Fab 평균 이물 개수는 `mean_particle_count / inspection_area_per_sample_mm2`로 밀도로 바꾸고,
  크기별 시뮬레이션 뒤 `exact_point`로 적용한다.

## 기하 규칙

- 이물은 직경 `d`, 반경 `r=d/2`인 전도성 구로 모델링한다. 중심 높이는
  `surface_z + r - embedding_depth`다.
- landing 시점까지 형성된 layer 접촉 여유는
  `r + contact_margin - hypot(Dxy, Dz)`다. 0 이상인 위치가 접촉 가능 위치다.
- 이후 형성 layer는 manifest의 `exclude`, `protrusion`, `wrap_contact` 정책을 따른다.
- `lateral_short`는 같은 layer의 다른 net, `vertical_short`는 다른 layer의 다른 net 조합이다.
  vertical 결과는 보정 전까지 `potential_vertical_short`다.
- 같은 조건의 전체 CA는 hotspot 지도를 OR한 뒤 적분한다. 전극 조합별 CA를 단순 합산하지 않는다.
- 겹친 픽셀의 설명용 기여도는 활성 불량 단위 수 `k`로 `1/k` 배분한다. short 기여는 두 전극에
  절반씩, open은 해당 전극에 전부 배분한다. raw와 unique CA도 함께 보존한다.

## 실행과 출력

- `simulate-nominal`: 명목 hotspot, CA, union, 내부 legacy kernel 참고값.
- `compare-sim`: 외부 기존 sim의 2-D binary 결과를 명시적 좌표 계약으로 읽어 공통 항목의 XOR와
  CA 차이를 계산한다. 비교 계약은 제품·설계 버전·mask hash가 현재 실행과 같아야 한다. 기존 sim에
  없던 전극 조합은 `new_pair`, open·CNO 단위는 `new_unit`, 적용 공정이 아닌 조건은 `not_applicable`이다.
- `simulate-probability`: Sobol scramble seed별 공정 산포 조건을 생성하고 조건별 판정 후 평균과
  p05/p50/p95를 계산한다. 여러 seed와 수렴 결과를 manifest에 남긴다.
- `apply-fab --run`: `risk_by_layer`, `risk_by_pair`, `risk_by_open_path`, `risk_by_electrode`,
  `risk_summary`, `risk_manifest`를 만든다.
- 모든 실행은 입력 hash, git commit/dirty, build ID, 단위, 가정, 결과 artifact hash를 남긴다.
- 대형 산포 실행은 예상 mask 메모리가 256 MB를 넘으면 packbits 임시 저장소를 사용한다.

## 개발 규칙과 남은 범위

Python 3.11+, CPU 전용이다. 입력 계약 위반은 `InputContractError`로 fail closed 처리한다. 사용자
메시지는 한국어, identifier와 파일 열 이름은 영어다. `ruff check`와 `pytest`를 통과해야 한다.

아직 물리 검증이 필요한 항목은 TR의 신호 분리 효과, 전도성/비전도성 이물별 contact margin,
원형 외 이물 형상·방향, vertical short 전기 보정, 비정규·상관 산포, Fab 원본 표본, tile+halo다.
