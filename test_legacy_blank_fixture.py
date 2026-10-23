"""历史 R50 叶片几何测试数据；不用于当前 Bezier 毛坯建模。"""

from dataclasses import dataclass


@dataclass(frozen=True)
class LegacyBlankParameters:
    hub_inlet_radius: float = 10.0
    front_inlet_radius: float = 27.5
    inlet_z: float = 40.0
    hub_bend_radius: float = 15.0
    front_bend_radius: float = 8.0
    hub_outlet_z: float = 0.0
    front_outlet_z: float = 10.0
    outlet_radius: float = 50.0

    @property
    def hub_center(self):
        return self.hub_inlet_radius+self.hub_bend_radius, self.hub_outlet_z+self.hub_bend_radius

    @property
    def front_center(self):
        return self.front_inlet_radius+self.front_bend_radius, self.front_outlet_z+self.front_bend_radius

    def validate(self):
        return self
