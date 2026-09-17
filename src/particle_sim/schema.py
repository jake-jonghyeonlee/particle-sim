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
    fab_distribution_layer: str


class OpenPath(BaseModel):
    """Open 판정 대상 전극과 단자. 단자가 없으면 fallback 규칙을 명시해야 한다."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    path_id: str
    electrode: str
    terminal_a_yx_px: tuple[int, int] | None = None
    terminal_b_yx_px: tuple[int, int] | None = None
    fallback: Literal["error", "roi_boundary", "skeleton_cut"] = "error"

    @model_validator(mode="after")
    def _check_terminals(self) -> OpenPath:
        given = (self.terminal_a_yx_px is None) + (self.terminal_b_yx_px is None)
        if given == 1:
            raise ValueError(
                f"open_path '{self.path_id}'의 단자는 둘 다 있거나 둘 다 없어야 합니다"
            )
        if given == 2 and self.fallback == "error":
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

    schema_version: Literal["1.0"]
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
    open_paths: list[OpenPath] = Field(default_factory=list)
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

        self._check_unique("open_paths", "path_id", [o.path_id for o in self.open_paths])
        for open_path in self.open_paths:
            if open_path.electrode not in electrode_by_id:
                raise ValueError(
                    f"open_path '{open_path.path_id}'의 전극 "
                    f"'{open_path.electrode}'이(가) 존재하지 않습니다"
                )

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

    schema_version: Literal["1.0"]
    assumption_set: str
    evidence: Literal["measured", "assumed"]
    thickness_delta_um: dict[str, NormalSpec] = Field(default_factory=dict)
    overlay_um: dict[str, NormalSpec] = Field(default_factory=dict)
    cd_bias_um: dict[str, NormalSpec] = Field(default_factory=dict)
    particle_diameter_delta_um: NormalSpec | None = None
    center_offset_um: NormalSpec | None = None
    embedding_depth_um: UniformSpec | None = None
    planarization_factor: dict[str, UniformSpec] = Field(default_factory=dict)
    contact_margin_um: UniformSpec | None = None

    def is_degenerate(self) -> bool:
        """모든 변수의 산포가 0이면 True. nominal 결과와 일치해야 한다."""
        normals = [
            *self.thickness_delta_um.values(),
            *self.overlay_um.values(),
            *self.cd_bias_um.values(),
        ]
        for spec in (self.particle_diameter_delta_um, self.center_offset_um):
            if spec is not None:
                normals.append(spec)
        if any(spec.sigma > 0.0 or spec.mean != 0.0 for spec in normals):
            return False
        uniforms = list(self.planarization_factor.values())
        for spec in (self.embedding_depth_um, self.contact_margin_um):
            if spec is not None:
                uniforms.append(spec)
        return all(spec.low == spec.high for spec in uniforms)


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
