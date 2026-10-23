"""双向差集体积核对；兼容服务端 IronPython 语法，不依赖 CAD 模块。"""


def overlap_from_losses(volume_a, volume_b, remaining_a, remaining_b, tolerance):
    """Compute intersection from symmetric difference losses; require agreement."""
    values = (volume_a, volume_b, remaining_a, remaining_b, tolerance)
    if any(value != value or abs(value) == float("inf") for value in values):
        raise ValueError("Non-finite volume in blade interference check")
    if volume_a <= 0 or volume_b <= 0 or tolerance <= 0:
        raise ValueError("Blade volumes and tolerance must be positive")
    if not (-tolerance <= remaining_a <= volume_a+tolerance
            and -tolerance <= remaining_b <= volume_b+tolerance):
        raise ValueError("Subtraction produced an invalid remaining volume")
    loss_a = volume_a-remaining_a
    loss_b = volume_b-remaining_b
    if abs(loss_a-loss_b) > tolerance:
        raise ValueError("Inconsistent overlap volume: losses %r, %r" % (loss_a, loss_b))
    overlap = max(0.0, (loss_a+loss_b)/2)
    if overlap > min(volume_a, volume_b)+tolerance:
        raise ValueError("Overlap exceeds the smaller blade volume")
    return overlap
