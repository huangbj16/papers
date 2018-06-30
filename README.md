# 光线追踪大作业
## 计63 黄冰鉴 2016011296
## 2018年6月

## 结果图
### Ray tracing
![ray_tracing](https://github.com/huangbj16/papers/blob/master/draw.png)
### Progressive photon mapping
![ppm](https://github.com/huangbj16/papers/blob/master/11.png)
![ppm](https://github.com/huangbj16/papers/blob/master/4.png)

## 实现内容
### 基本效果：反射，折射，阴影。
```C++
//反射
if (thing[p]->material.reflection > 0) {
  Line newline;
  newline.dir = line.dir.Reflect(vecN);
  newline.point = crash_bag.crash_point + newline.dir * EPS;
  colorlight += Intersect(newline, time + 1, i, j, weight * thing[p]->material.reflection) * thing[p]->material.Getcolor(crash_bag.u, crash_bag.v) * thing[p]->material.reflection;//calculate reflect color
}
//折射
if (thing[p]->material.refraction > 0) {
  double n = thing[p]->material.refraction_index;
  Line newline;
  n = 1 / n;
  newline.dir = line.dir.Refract(vecN, n);
  if (!newline.dir.IsZeroVector()) {
    newline.point = crash_bag.crash_point + newline.dir * EPS;
    double absorb = exp((-0.01f) * thing[p]->LengthInside(vecN, newline.dir));
    colorlight += Intersect(newline, time + 1, i, j, weight * thing[p]->material.refraction * absorb) * (thing[p]->material.refraction * absorb);
  }
}
//阴影
for (int i = 0; i < thingnum; ++i) {
  if (i == p) continue;
  other_bag = thing[i]->Crash(crash_point, l);
  if (!other_bag.crash_point.IsNullVector()) {
    if ((other_bag.crash_point - crash_point).Module2() < dist) {
      return Color(0, 0, 0);
    }
  }
}
//phong模型漫反射和高光
if (crash_thing->material.diffusion > 0) {
  double dot = l.Dot(vecN);
  if (dot > 0) {
    double diff = dot * crash_thing->material.diffusion;
    c += crash_thing->material.Getcolor(bag.u, bag.v) * crash_light->material.color * diff;
  }
}
if (crash_thing->material.specular > 0) {
  double dot = r.Dot(view_direction);
  if (dot > 0) {
    c += crash_light->material.color * crash_thing->material.specular * pow(dot, 20);
  }
}
```




结果图，实现的内容及对应代码段
:
简述实现过程，比如牛顿迭代
