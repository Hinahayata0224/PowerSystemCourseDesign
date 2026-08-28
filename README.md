# Power System Course Design (C2)

基于 MATLAB 和 DSP（PSD-BPA）的**输电网稳定性计算分析**课程设计。

## 目录结构

```
Power System/
├── 报告与演示/     # 最终交付：设计报告、汇报 PPT、课程设计文档
├── 源代码/         # 最终版 MATLAB 脚本与 PSD-BPA 数据/故障文件
└── archive/        # 过程文件归档：个人报告、故障说明、报告模板、DSP 工程等
```

## 内容概览

- **任务**：输电网潮流计算、短路计算、稳定性分析（C2 课题）
- **工具**：MATLAB（潮流/短路脚本）+ PSD-BPA（DSP 稳定计算）
- **成员**：林赵鑫（202230232118）、许啸扬（202230234075）

## 源代码清单

| 文件 | 说明 |
|---|---|
| `C2task1_dailyload.m` | 日负荷曲线 |
| `C2task2_powerflow.m` | 潮流计算（主） |
| `C2task2_powerflow_single56.m` | 潮流计算（单线 5-6） |
| `C2task3_powerflow.dat` | 潮流数据 |
| `C2task56_shortcircuit_node6_ode45.m` | 节点 6 短路（ode45） |
| `C2task7_*.swi` | PSD-BPA 故障/短路稳定文件 |

> archive/ 为过程版本归档，非最终交付。
