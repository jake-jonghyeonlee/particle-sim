"""pydantic v2 기반 ProductManifest 스키마."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

MaterialType = Literal["conductor", "inorganic", "organic", "other"]
SurfaceModel = Literal["flat", "conformal", "partial"]
FailureMode = Literal["lateral_short", "vertical_short"]

PARTICLE_SIZES_UM: tuple[float, ...] = tuple(float(i) for i in range(1, 31)) + (40.0, 50.0)


class CoordinateSystem(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    origin: Literal["upper_left"]
    x_direction: Literal["right"]
    y_direction: Literal["down"]
    pixel_size_x_um: float = Field(gt=0)
    pixel_size_y_um: float = Field(gt=0)
    roi_width_um: float = Field(gt=0)
    roi_height_um: float = Field(gt=0)


class StackLayer(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    layer_name: str
    process_order: int = Field(ge=0)
    material_type: MaterialType
    thickness_um: float = Field(gt=0)
    thickness_sigma_um: float | None = Field(default=None, ge=0)
    planarization_factor: float | None = Field(default=None, ge=0, le=1)
    z_bottom_um: float | None = None
    patterned: bool


class LayerMask(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    layer_name: str
    path: str


class Electrode(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    electrode_id: str
    layer_name: str
    net_id: str
    electrode_role: str
    path: str


class ParticleSourceStep(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    step_id: str
    lands_after_layer: str
    surface_model: SurfaceModel
    # 생략하면 lands_after_layer 바로 위 마스크 layer(patterned 또는 contact hole layer) 이름으로
    # 정한다. schema 2.0에서 명시하면 그 이름과 같아야 하고, 위에 마스크 layer가 없을 때만 적는다.
    fab_distribution_layer: str | None = None


class OpenPath(BaseModel):
    """Open 판정 대상 전극과 단자. 단자가 없으면 fallback 규칙을 명시해야 한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    path_id: str
    electrode: str
    terminal_a_yx_px: tuple[int, int] | None = None
    terminal_b_yx_px: tuple[int, int] | None = None
    terminal_a_mask_path: str | None = None
    terminal_b_mask_path: str | None = None
    applicable_particle_source_steps: list[str] = Field(default_factory=list)
    # skeleton_cut(중심선 절단 proxy)은 아직 구현하지 않아 스키마에서도 허용하지 않는다
    fallback: Literal["error", "roi_boundary"] = "error"

    @model_validator(mode="after")
    def _check_terminals(self) -> OpenPath:
        point_count = int(self.terminal_a_yx_px is not None) + int(
            self.terminal_b_yx_px is not None
        )
        mask_count = int(self.terminal_a_mask_path is not None) + int(
            self.terminal_b_mask_path is not None
        )
        if point_count == 1:
            raise ValueError(
                f"open_path '{self.path_id}'의 단자는 둘 다 있거나 둘 다 없어야 합니다"
            )
        if mask_count == 1:
            raise ValueError(
                f"open_path '{self.path_id}'의 단자 mask는 둘 다 있거나 둘 다 없어야 합니다"
            )
        if point_count == 2 and mask_count == 2:
            raise ValueError(
                f"open_path '{self.path_id}'는 단자 좌표와 단자 mask를 함께 쓸 수 없습니다"
            )
        if point_count == 0 and mask_count == 0 and self.fallback == "error":
            raise ValueError(
                f"open_path '{self.path_id}'에 단자가 없으면 fallback을 지정해야 합니다"
            )
        for name, terminal in (
            ("terminal_a_yx_px", self.terminal_a_yx_px),
            ("terminal_b_yx_px", self.terminal_b_yx_px),
        ):
            if terminal is not None and (terminal[0] < 0 or terminal[1] < 0):
                raise ValueError(f"open_path '{self.path_id}'의 {name}는 0 이상이어야 합니다")
        return self


AUTO_OPEN_PREFIX = "auto_open:"
CNO_PREFIX = "cno:"
RESERVED_UNIT_PREFIXES = (AUTO_OPEN_PREFIX, CNO_PREFIX)


class Contact(BaseModel):
    """contact·via 원본 mask. 각 연결 성분이 위·아래 도체 전극을 잇는 hole 하나다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    contact_id: str = Field(min_length=1)
    path: str
    hole_layer: str
    # 생략하면 stack process_order에서 hole_layer 바로 아래·위의 patterned 도체 layer를 쓴다
    lower_layer: str | None = None
    upper_layer: str | None = None
    applicable_particle_source_steps: list[str] = Field(default_factory=list)


def resolve_contact_layers(
    stack: list[StackLayer], contact: Contact
) -> tuple[str | None, str | None]:
    """contact가 잇는 아래·위 도체 layer. 찾지 못하면 None."""
    ordered = sorted(stack, key=lambda layer: layer.process_order)
    names = [layer.layer_name for layer in ordered]
    if contact.hole_layer not in names:
        return contact.lower_layer, contact.upper_layer
    hole_index = names.index(contact.hole_layer)

    def is_conductor(layer: StackLayer) -> bool:
        return layer.patterned and layer.material_type == "conductor"

    lower = contact.lower_layer
    if lower is None:
        below = [layer for layer in ordered[:hole_index] if is_conductor(layer)]
        lower = below[-1].layer_name if below else None
    upper = contact.upper_layer
    if upper is None:
        above = [layer for layer in ordered[hole_index + 1 :] if is_conductor(layer)]
        upper = above[0].layer_name if above else None
    return lower, upper


class ExcludedElectrodePair(BaseModel):
    """계산에서 제외하기로 승인한 다른 net 전극 조합."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    electrode_a: str
    electrode_b: str
    reason: str = Field(min_length=1)

    @field_validator("reason")
    @classmethod
    def _reason_not_blank(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("제외 전극 조합의 reason은 공백일 수 없습니다")
        return value


class ElectrodePair(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    pair_id: str
    electrode_a: str
    electrode_b: str
    failure_mode: FailureMode


class LegacyKernel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    size: int = Field(default=5, ge=3)
    shape: Literal["square", "disk"] = "square"
    iterations: Literal["round", "ceil", "floor"] = "round"

    @model_validator(mode="after")
    def _check_odd_size(self) -> LegacyKernel:
        if self.size % 2 == 0:
            raise ValueError(f"legacy_kernel.size는 3 이상의 홀수여야 합니다: {self.size}")
        return self


class PostStepPolicy(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    mode: Literal["exclude", "protrusion", "wrap_contact"]
    coverage_break_factor: float | None = Field(default=None, gt=0)

    @model_validator(mode="after")
    def _check_coverage_break_factor(self) -> PostStepPolicy:
        if self.mode == "protrusion" and self.coverage_break_factor is None:
            raise ValueError("protrusion 모드에는 coverage_break_factor가 필요합니다")
        return self


class ProductManifest(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: Literal["1.0", "2.0"]
    product_id: str
    product_family: str
    design_version: str
    design_lineage_id: str
    fab_id: str
    is_synthetic: bool
    coordinate_system: CoordinateSystem
    stack: list[StackLayer]
    layer_masks: list[LayerMask]
    electrodes: list[Electrode]
    particle_source_steps: list[ParticleSourceStep]
    electrode_pairs: list[ElectrodePair]
    excluded_electrode_pairs: list[ExcludedElectrodePair] = Field(default_factory=list)
    open_paths: list[OpenPath] = Field(default_factory=list)
    contacts: list[Contact] = Field(default_factory=list)
    # contact·ROI 경계 자동 단자 open을 계산할 layer와 적용 이물 유입 공정
    open_applicable_steps_by_layer: dict[str, list[str]] = Field(default_factory=dict)
    particle_sizes_um: list[float] = Field(default_factory=lambda: list(PARTICLE_SIZES_UM))
    embedding_depth_um: float = Field(ge=0)
    contact_margin_um: float
    legacy_kernel: LegacyKernel = Field(default_factory=LegacyKernel)
    post_step_policy: PostStepPolicy

    @field_validator("particle_sizes_um")
    @classmethod
    def _check_particle_sizes(cls, value: list[float]) -> list[float]:
        if not value:
            raise ValueError("particle_sizes_um은 비어 있을 수 없습니다")
        if any(size <= 0 for size in value):
            raise ValueError("particle_sizes_um은 모두 0보다 커야 합니다")
        if any(b <= a for a, b in zip(value, value[1:], strict=False)):
            raise ValueError("particle_sizes_um은 엄격히 증가해야 합니다")
        return value

    def upper_mask_layer(self, layer_name: str) -> str | None:
        """process_order에서 layer_name 바로 위의 마스크 layer. 없으면 None.

        마스크 layer는 patterned layer와 contact hole layer(contact mask가 있는 절연막)다.
        마스크 없는 절연막은 건너뛴다.
        """
        hole_layers = {contact.hole_layer for contact in self.contacts}
        ordered = sorted(self.stack, key=lambda layer: layer.process_order)
        order = next(x.process_order for x in ordered if x.layer_name == layer_name)
        above = [
            x for x in ordered
            if x.process_order > order and (x.patterned or x.layer_name in hole_layers)
        ]
        return above[0].layer_name if above else None

    def _resolve_fab_layer(self, step: ParticleSourceStep) -> ParticleSourceStep:
        """Fab 이물 밀도는 착지 layer 바로 위 마스크 layer의 값을 쓴다(2026-09-30 결정)."""
        upper = self.upper_mask_layer(step.lands_after_layer)
        given = step.fab_distribution_layer
        if given and self.schema_version != "2.0":
            return step  # schema 1.0의 명시값은 이전처럼 자유 식별자로 둔다
        if upper is None:
            if not given:
                raise ValueError(
                    f"step '{step.step_id}': lands_after_layer '{step.lands_after_layer}' 위에 "
                    "마스크 layer가 없어 fab_distribution_layer를 명시해야 합니다"
                )
            return step
        if given and given != upper:
            raise ValueError(
                f"step '{step.step_id}': fab_distribution_layer '{given}'은(는) "
                f"lands_after_layer '{step.lands_after_layer}' 바로 위 마스크 layer "
                f"'{upper}'와 같아야 합니다"
            )
        return step if given else step.model_copy(update={"fab_distribution_layer": upper})

    @model_validator(mode="after")
    def _check_cross_references(self) -> ProductManifest:
        stack_by_name = {layer.layer_name: layer for layer in self.stack}
        self._check_unique("stack", "layer_name", [x.layer_name for x in self.stack])
        self._check_unique("stack", "process_order", [x.process_order for x in self.stack])
        self._check_unique("electrodes", "electrode_id", [e.electrode_id for e in self.electrodes])
        self._check_unique(
            "particle_source_steps", "step_id", [s.step_id for s in self.particle_source_steps]
        )
        self._check_unique("electrode_pairs", "pair_id", [p.pair_id for p in self.electrode_pairs])

        layer_names = set(stack_by_name)
        for lm in self.layer_masks:
            if lm.layer_name not in layer_names:
                raise ValueError(
                    f"layer_masks의 layer_name '{lm.layer_name}'이(가) stack에 없습니다"
                )
        for el in self.electrodes:
            if el.layer_name not in layer_names:
                raise ValueError(
                    f"electrodes의 layer_name '{el.layer_name}'이(가) stack에 없습니다"
                )

        patterned = {x.layer_name for x in self.stack if x.patterned}
        mask_names = [lm.layer_name for lm in self.layer_masks]
        for name in sorted(patterned):
            if mask_names.count(name) != 1:
                raise ValueError(
                    f"patterned 레이어 '{name}'은(는) layer_masks에 정확히 하나 있어야 합니다"
                )
        for name in sorted(layer_names - patterned):
            if name in mask_names:
                raise ValueError(f"non-patterned 레이어 '{name}'은(는) layer_masks에 없어야 합니다")

        for el in self.electrodes:
            if el.layer_name not in patterned:
                raise ValueError(
                    f"electrode '{el.electrode_id}'의 레이어 '{el.layer_name}'은(는) "
                    "patterned여야 합니다"
                )
            if stack_by_name[el.layer_name].material_type != "conductor":
                raise ValueError(
                    f"electrode '{el.electrode_id}'의 레이어 '{el.layer_name}'은(는) "
                    "conductor여야 합니다"
                )

        for step in self.particle_source_steps:
            if step.lands_after_layer not in layer_names:
                raise ValueError(
                    f"step '{step.step_id}'의 lands_after_layer "
                    f"'{step.lands_after_layer}'이(가) stack에 없습니다"
                )
        object.__setattr__(
            self,
            "particle_source_steps",
            [self._resolve_fab_layer(step) for step in self.particle_source_steps],
        )

        electrode_by_id = {e.electrode_id: e for e in self.electrodes}
        for pair in self.electrode_pairs:
            for side, eid in (("electrode_a", pair.electrode_a), ("electrode_b", pair.electrode_b)):
                if eid not in electrode_by_id:
                    raise ValueError(
                        f"pair '{pair.pair_id}'의 {side} '{eid}'이(가) 존재하지 않습니다"
                    )
            if pair.electrode_a == pair.electrode_b:
                raise ValueError(
                    f"pair '{pair.pair_id}'의 electrode_a와 electrode_b는 달라야 합니다"
                )
            ea = electrode_by_id[pair.electrode_a]
            eb = electrode_by_id[pair.electrode_b]
            if ea.net_id == eb.net_id:
                raise ValueError(f"pair '{pair.pair_id}'의 net_id는 달라야 합니다")
            if pair.failure_mode == "lateral_short" and ea.layer_name != eb.layer_name:
                raise ValueError(
                    f"lateral_short pair '{pair.pair_id}'은(는) 같은 layer_name이어야 합니다"
                )
            if pair.failure_mode == "vertical_short" and ea.layer_name == eb.layer_name:
                raise ValueError(
                    f"vertical_short pair '{pair.pair_id}'은(는) 다른 layer_name이어야 합니다"
                )

        declared_keys: set[tuple[str, str]] = set()
        for pair in self.electrode_pairs:
            key = tuple(sorted((pair.electrode_a, pair.electrode_b)))
            if key in declared_keys:
                raise ValueError(f"전극 조합 {key}이(가) 다른 pair_id로 중복 선언됐습니다")
            declared_keys.add(key)
        excluded_keys: set[tuple[str, str]] = set()
        for excluded in self.excluded_electrode_pairs:
            if excluded.electrode_a not in electrode_by_id:
                raise ValueError(
                    f"제외 조합의 electrode_a '{excluded.electrode_a}'이(가) 존재하지 않습니다"
                )
            if excluded.electrode_b not in electrode_by_id:
                raise ValueError(
                    f"제외 조합의 electrode_b '{excluded.electrode_b}'이(가) 존재하지 않습니다"
                )
            if excluded.electrode_a == excluded.electrode_b:
                raise ValueError("제외 조합의 두 전극은 달라야 합니다")
            ea = electrode_by_id[excluded.electrode_a]
            eb = electrode_by_id[excluded.electrode_b]
            if ea.net_id == eb.net_id:
                raise ValueError("같은 net 전극은 제외 조합에 넣지 않습니다")
            key = tuple(sorted((excluded.electrode_a, excluded.electrode_b)))
            if key in excluded_keys:
                raise ValueError(f"제외 전극 조합 {key}이(가) 중복됩니다")
            if key in declared_keys:
                raise ValueError(f"전극 조합 {key}이(가) 계산과 제외에 동시에 선언됐습니다")
            excluded_keys.add(key)

        self._check_unique("open_paths", "path_id", [o.path_id for o in self.open_paths])
        # 결과 배열과 배분표는 불량 단위 ID로 묶이므로 pair와 open path가 ID를 공유하면 안 된다.
        shared_ids = sorted(
            {p.pair_id for p in self.electrode_pairs} & {o.path_id for o in self.open_paths}
        )
        if shared_ids:
            raise ValueError(
                f"electrode_pairs의 pair_id와 open_paths의 path_id가 겹칩니다: {shared_ids}"
            )
        for open_path in self.open_paths:
            if open_path.electrode not in electrode_by_id:
                raise ValueError(
                    f"open_path '{open_path.path_id}'의 전극 "
                    f"'{open_path.electrode}'이(가) 존재하지 않습니다"
                )
            unknown_steps = sorted(
                set(open_path.applicable_particle_source_steps)
                - {step.step_id for step in self.particle_source_steps}
            )
            if unknown_steps:
                raise ValueError(
                    f"open_path '{open_path.path_id}'의 적용 이물 유입 공정이 존재하지 않습니다: "
                    f"{unknown_steps}"
                )
            if self.schema_version == "2.0" and not open_path.applicable_particle_source_steps:
                raise ValueError(
                    f"schema 2.0 open_path '{open_path.path_id}'에는 "
                    "applicable_particle_source_steps가 필요합니다"
                )
            if len(open_path.applicable_particle_source_steps) != len(
                set(open_path.applicable_particle_source_steps)
            ):
                raise ValueError(
                    f"open_path '{open_path.path_id}'의 적용 이물 유입 공정이 중복됩니다"
                )

        self._check_contacts(stack_by_name)

        for step in self.particle_source_steps:
            if step.surface_model != "partial":
                continue
            order = stack_by_name[step.lands_after_layer].process_order
            for layer in self.stack:
                if (
                    layer.process_order <= order
                    and layer.material_type != "conductor"
                    and layer.planarization_factor is None
                ):
                    raise ValueError(
                        f"partial surface_model step '{step.step_id}'에서 레이어 "
                        f"'{layer.layer_name}'은(는) planarization_factor가 필요합니다"
                    )

        if self.post_step_policy.mode == "protrusion":
            for layer in self.stack:
                if layer.material_type != "conductor" and layer.planarization_factor is None:
                    raise ValueError(
                        f"protrusion 모드에서 레이어 '{layer.layer_name}'은(는) "
                        "planarization_factor가 필요합니다"
                    )

        return self

    def _check_contacts(self, stack_by_name: dict[str, StackLayer]) -> None:
        unit_ids = [p.pair_id for p in self.electrode_pairs] + [
            o.path_id for o in self.open_paths
        ]
        for unit_id in unit_ids:
            if unit_id.startswith(RESERVED_UNIT_PREFIXES):
                raise ValueError(
                    f"불량 단위 ID '{unit_id}'은(는) 예약 접두사 {RESERVED_UNIT_PREFIXES}를 "
                    "쓸 수 없습니다"
                )
        if (self.contacts or self.open_applicable_steps_by_layer) and self.schema_version != "2.0":
            raise ValueError("contacts와 open_applicable_steps_by_layer는 schema 2.0에서만 씁니다")
        step_ids = {step.step_id for step in self.particle_source_steps}
        electrode_layers = {e.layer_name for e in self.electrodes}
        self._check_unique("contacts", "contact_id", [c.contact_id for c in self.contacts])
        contact_layers: set[str] = set()
        for contact in self.contacts:
            name = f"contact '{contact.contact_id}'"
            if ":" in contact.contact_id or "__" in contact.contact_id:
                raise ValueError(f"{name}의 contact_id에는 ':'와 '__'를 쓸 수 없습니다")
            hole = stack_by_name.get(contact.hole_layer)
            if hole is None:
                raise ValueError(f"{name}의 hole_layer '{contact.hole_layer}'이(가) stack에 없습니다")
            if hole.material_type == "conductor":
                raise ValueError(f"{name}의 hole_layer는 절연층이어야 합니다")
            lower, upper = resolve_contact_layers(self.stack, contact)
            for side, layer_name in (("lower_layer", lower), ("upper_layer", upper)):
                if layer_name is None:
                    raise ValueError(f"{name}의 {side}를 stack에서 찾을 수 없습니다")
                layer = stack_by_name.get(layer_name)
                if layer is None or not (layer.patterned and layer.material_type == "conductor"):
                    raise ValueError(
                        f"{name}의 {side} '{layer_name}'은(는) patterned 도체 layer여야 합니다"
                    )
                if layer_name not in electrode_layers:
                    raise ValueError(f"{name}의 {side} '{layer_name}'에 전극이 없습니다")
            assert lower is not None and upper is not None
            order = stack_by_name[contact.hole_layer].process_order
            if not stack_by_name[lower].process_order < order < stack_by_name[upper].process_order:
                raise ValueError(f"{name}의 hole_layer는 lower_layer와 upper_layer 사이여야 합니다")
            if not contact.applicable_particle_source_steps:
                raise ValueError(f"{name}에는 applicable_particle_source_steps가 필요합니다")
            unknown = sorted(set(contact.applicable_particle_source_steps) - step_ids)
            if unknown:
                raise ValueError(f"{name}의 적용 이물 유입 공정이 존재하지 않습니다: {unknown}")
            if len(set(contact.applicable_particle_source_steps)) != len(
                contact.applicable_particle_source_steps
            ):
                raise ValueError(f"{name}의 적용 이물 유입 공정이 중복됩니다")
            contact_layers.update((lower, upper))
        for layer_name, steps in self.open_applicable_steps_by_layer.items():
            if layer_name not in electrode_layers:
                raise ValueError(
                    f"open_applicable_steps_by_layer의 '{layer_name}'은(는) 전극이 있는 layer여야 합니다"
                )
            unknown = sorted(set(steps) - step_ids)
            if unknown:
                raise ValueError(
                    f"open_applicable_steps_by_layer['{layer_name}']의 공정이 존재하지 않습니다: "
                    f"{unknown}"
                )
            if len(set(steps)) != len(steps):
                raise ValueError(f"open_applicable_steps_by_layer['{layer_name}']가 중복됩니다")
        missing = sorted(contact_layers - set(self.open_applicable_steps_by_layer))
        if missing:
            raise ValueError(
                "contact로 이어진 layer는 open_applicable_steps_by_layer에 선언해야 합니다"
                f"(자동 open을 계산하지 않으면 빈 목록): {missing}"
            )

    @staticmethod
    def _check_unique(section: str, key: str, values: list[object]) -> None:
        seen: set[object] = set()
        for v in values:
            if v in seen:
                raise ValueError(f"{section}의 {key} '{v}'이(가) 중복됩니다")
            seen.add(v)


class NormalSpec(BaseModel):
    """평균과 표준편차로 정의한 정규 분포."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    mean: float = 0.0
    sigma: float = Field(ge=0.0)
    evidence: Literal["measured", "assumed"] | None = None
    source_id: str | None = None
    source_hash: str | None = None
    measurement_period: str | None = None
    fab_process: str | None = None
    product_family: str | None = None


class OverlayXYSpec(BaseModel):
    """x/y가 서로 다른 overlay 분포."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    x: NormalSpec
    y: NormalSpec


class UniformSpec(BaseModel):
    """구간 균등 분포. low == high면 고정값이다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    low: float
    high: float

    @model_validator(mode="after")
    def _check_range(self) -> UniformSpec:
        if self.high < self.low:
            raise ValueError(f"uniform 구간이 뒤집혔습니다: low={self.low}, high={self.high}")
        return self


class ScenarioConfig(BaseModel):
    """공정 산포 시나리오 설정. 실제 측정값이 아니면 evidence='assumed'로 기록한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: Literal["1.0", "2.0"]
    assumption_set: str
    evidence: Literal["measured", "assumed"]
    thickness_delta_um: dict[str, NormalSpec] = Field(default_factory=dict)
    overlay_um: dict[str, NormalSpec | OverlayXYSpec] = Field(default_factory=dict)
    cd_bias_um: dict[str, NormalSpec] = Field(default_factory=dict)
    overlay_reference_layer: str | None = None
    overlay_sign_convention: Literal["positive_x_right_positive_y_down"] | None = None
    independence_model: Literal["independent"] = "independent"
    particle_diameter_delta_um: NormalSpec | None = None
    center_offset_um: NormalSpec | None = None
    embedding_depth_um: UniformSpec | None = None
    planarization_factor: dict[str, UniformSpec] = Field(default_factory=dict)
    contact_margin_um: UniformSpec | None = None

    def is_degenerate(self) -> bool:
        """산포가 모두 0이면 True(평균이 0이 아닌 고정 이동은 산포가 아니므로 허용한다)."""
        overlay_normals = []
        for spec in self.overlay_um.values():
            if isinstance(spec, OverlayXYSpec):
                overlay_normals.extend((spec.x, spec.y))
            else:
                overlay_normals.extend((spec, spec))
        normals = [*self.thickness_delta_um.values(), *overlay_normals, *self.cd_bias_um.values()]
        for spec in (self.particle_diameter_delta_um, self.center_offset_um):
            if spec is not None:
                normals.append(spec)
        if any(spec.sigma > 0.0 for spec in normals):
            return False
        uniforms = list(self.planarization_factor.values())
        for spec in (self.embedding_depth_um, self.contact_margin_um):
            if spec is not None:
                uniforms.append(spec)
        return all(spec.low == spec.high for spec in uniforms)

    @model_validator(mode="after")
    def _check_v2_overlay_contract(self) -> ScenarioConfig:
        if self.schema_version != "2.0":
            return self
        legacy_layers = [
            layer for layer, spec in self.overlay_um.items() if isinstance(spec, NormalSpec)
        ]
        if legacy_layers:
            raise ValueError(
                f"schema 2.0 overlay는 x/y 분포를 각각 지정해야 합니다: {legacy_layers}"
            )
        if self.overlay_um and self.overlay_reference_layer is None:
            raise ValueError("schema 2.0 overlay에는 overlay_reference_layer가 필요합니다")
        if self.overlay_um and self.overlay_sign_convention is None:
            raise ValueError("schema 2.0 overlay에는 overlay_sign_convention이 필요합니다")
        measured_specs: list[tuple[str, NormalSpec]] = []
        for layer, spec in self.overlay_um.items():
            assert isinstance(spec, OverlayXYSpec)
            measured_specs.extend(((f"overlay_um.{layer}.x", spec.x), (f"overlay_um.{layer}.y", spec.y)))
        measured_specs.extend(
            (f"cd_bias_um.{layer}", spec) for layer, spec in self.cd_bias_um.items()
        )
        for name, spec in measured_specs:
            if spec.evidence is None:
                raise ValueError(f"schema 2.0 {name}에는 evidence가 필요합니다")
            if spec.evidence == "measured" and spec.source_id is None:
                raise ValueError(f"실측 {name}에는 source_id가 필요합니다")
            if spec.source_id is not None and not spec.source_id.strip():
                raise ValueError(f"{name}의 source_id는 공백일 수 없습니다")
        return self


def scenario_warnings(config: ScenarioConfig) -> list[str]:
    """시나리오 설정의 근거를 TODO_PHYSICS 경고로 만든다."""
    warnings = [
        f"TODO_PHYSICS: scenario.assumption_set={config.assumption_set} "
        f"evidence={config.evidence} 산포 입력의 근거입니다"
    ]
    if config.evidence == "assumed":
        warnings.append(
            "TODO_PHYSICS: scenario.evidence=assumed 확률값은 가정한 산포에 조건부이며 "
            "실제 Fab 불량확률이 아닙니다"
        )
    return warnings


def todo_physics_warnings(m: ProductManifest) -> list[str]:
    """TODO_PHYSICS 경고 목록을 생성한다."""
    policy = m.post_step_policy
    warnings: list[str] = [
        "TODO_PHYSICS: particle_conductivity=conductive 이물은 전도성으로 가정합니다",
        f"TODO_PHYSICS: embedding_depth_um={m.embedding_depth_um} 매립 깊이는 미확정 물리값입니다",
        f"TODO_PHYSICS: contact_margin_um={m.contact_margin_um} 접촉 여유값은 미확정 물리값입니다",
        f"TODO_PHYSICS: post_step_policy.mode={policy.mode} 이후 형성 layer 접촉 정책은 가정입니다",
        "TODO_PHYSICS: electrode_z=nominal 전극 높이는 하부 topography와 무관한 평균값입니다",
    ]
    if any(m.open_applicable_steps_by_layer.values()):
        warnings.append(
            "TODO_PHYSICS: 자동 단자 open은 전극 단위 판정이며 같은 net의 다른 layer 우회 경로는 "
            "반영하지 않습니다"
        )
    for open_path in m.open_paths:
        if not open_path.applicable_particle_source_steps:
            warnings.append(
                f"TODO_PHYSICS: open_path={open_path.path_id}의 적용 이물 유입 공정이 없어 "
                "v1 호환 동작으로 모든 공정에 계산합니다"
            )
    if policy.coverage_break_factor is not None:
        warnings.append(
            f"TODO_PHYSICS: post_step_policy.coverage_break_factor={policy.coverage_break_factor} "
            "절연막 피복 파괴 계수는 미확정 물리값입니다"
        )
    for layer in m.stack:
        if layer.planarization_factor is not None:
            warnings.append(
                f"TODO_PHYSICS: {layer.layer_name}.planarization_factor="
                f"{layer.planarization_factor} 평탄화 계수는 미확정 물리값입니다"
            )
    return warnings
