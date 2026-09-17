"""particle-sim 에러 타입."""

from __future__ import annotations


class InputContractError(ValueError):
    """입력 계약 위반 에러.

    str(e)는 "[product=P][file=F][field=X] message" 형식이며,
    값이 없는 항목은 생략된다.
    """

    def __init__(
        self,
        message: str,
        *,
        product_id: str | None = None,
        file: str | None = None,
        field: str | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.product_id = product_id
        self.file = file
        self.field = field

    def __str__(self) -> str:
        parts: list[str] = []
        if self.product_id is not None:
            parts.append(f"[product={self.product_id}]")
        if self.file is not None:
            parts.append(f"[file={self.file}]")
        if self.field is not None:
            parts.append(f"[field={self.field}]")
        prefix = "".join(parts)
        return f"{prefix} {self.message}" if prefix else self.message
