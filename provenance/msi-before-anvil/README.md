# IntegMultiReg → UMN MSI 运行与交接包

日期：2026-09-23。这个包用于在 Linux 集群上重新安装并运行研究任务，不是最终论文或 CRAN 发布包。

## 上传后先做什么

将整个 `.tar.gz` 和同名 `.sha256` 上传到 MSI 的持久项目目录。不要把长期输出放在 `/tmp`。建议预留至少 30 GB；完整模拟的任务原件和汇总副本会同时保留，实际占用以运行后为准。

```bash
sha256sum -c IntegMultiReg-MSI-20260923.tar.gz.sha256
tar -xzf IntegMultiReg-MSI-20260923.tar.gz
cd IntegMultiReg-MSI-20260923
sha256sum -c SHA256SUMS
cp msi/config.example.sh msi/config.sh
module spider R
module spider gsl
# 编辑 msi/config.sh：填写 account、实际 R/GSL/编译器模块名和持久输出路径
bash msi/submit.sh setup
```

使用 **R ≥ 4.4、GSL ≥ 2、C/C++/Fortran 编译器、make、unzip、Python 3、flock**。源码依赖已随包提供，安装不需要从 CRAN 下载。不要在 login node 上运行安装或采样；setup 会作为 Slurm 作业在计算节点上安装。记录它返回的 job ID，通过以下命令确认完成：

```bash
squeue --me
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS
# setup 成功后（library-linux/READY 存在）才提交 smoke
bash msi/submit.sh smoke
# smoke 成功后（runs/smoke/PASS 存在）才提交研究计算
bash msi/submit.sh chains
```

日志位于配置的 `MSI_RESULTS/logs/`。默认 `msismall`、每任务 1 CPU / 16 GB / 24 小时，数组最多同时运行 4 个任务。这些是起始资源申请，不能保证每条链都在限额内完成；用 `sacct` 核查真实耗时和峰值内存。account 和 module 必须使用你在 MSI 的实际配置，本包没有猜测它们。

## 任务和科学设置

| 命令 | 作业数量 | 内容 |
|---|---:|---|
| `chains` | 8 | 每条 1,000,000 retained draws + 50,000 burn-in；paper sampler、临床协变量、nu=(-4,-3,-4)、seed=100–107、八组显式起点 |
| `simulation` | 150 | 3 配置 × 50 次重复；每任务一个配置/重复，包含 IMR、BMS、L1-CPH、Uni-CPH |
| `correlated` | 180 | 6 配置 × 30 次重复，方法同上 |

模拟固定为 **code2017 / fixed markers**、350,000 retained draws、50,000 burn-in、10 folds、1 round、workers=1。第 r 次重复的 seed=100+r−1；不同配置共享该重复编号的 seed，与现有 runner 一致。只有决定执行这条已标明的方法路线时，才提交：

```bash
bash msi/submit.sh simulation
bash msi/submit.sh correlated
```

这不会解决 paper-reference 随机 marker / 无惩罚 CV 的秩亏问题。Training-refit 或 ridge sensitivity 需要另立协议；此包没有静默替换方法，也没有启动这些替代实验。

Linux 上的长链是新的、同设定的诊断运行。R/GSL/编译器和浮点环境变化可能使轨迹不同，**不能把它与 Mac 的前 350,000 draws 当作逐位相同的延长链**。保留旧八链作为对照。若需要精确前缀检验，应在同一个 MSI 安装环境另跑原长度基线，并使用 `paper/check-chain-prefix.R` 的严格检查，不能削弱其环境与 hash 条件。

## 完成后汇总和验收

确认对应数组任务都成功后，在 login node 提交汇总作业：

```bash
bash msi/submit.sh collect chains
bash msi/submit.sh collect simulation
bash msi/submit.sh collect correlated
```

汇总器要求所有任务的全局参数、源码与数据 hash 完全一致，并验证任务编号和完成状态；缺任务、参数混用或失败任务会被拒绝。原文件保存在 `runs/TYPE/tasks/`，汇总副本在 `runs/TYPE/combined/`。链汇总检查完整 checkpoint 与诊断投影一致性，并生成 PSRF、rank Rhat、ESS、MCSE、排名稳定性、前后半段漂移和诊断图。模拟汇总调用完整 50/30 次重复验收、baseline 验收和重复间 SD/SE 汇总。

`VALIDATION-COMPLETED` 只说明这些程序检查通过。仍需审阅警告、混合/精度和图表，才能更新论文结论；不能把程序完成等同于收敛或历史表格逐位复现。

## 失败、重跑和环境冻结

每个任务独立输出并持有进程锁；同一任务不能同时写入。重新提交相同类型的数组会跳过匹配设置且已完成的结果。现有 sampler **没有 native mid-chain checkpoint**：未完成的一条链需要从头重跑；已经写完的 fit 可以用于后续恢复。

正常退出会写入 `slurm-JOBID-TASK.exitcode`。节点故障、OOM、SIGKILL 可能来不及写退出文件，R 的 status 也可能停在 running；这时以 `sacct` 的终态及 Slurm 日志为准，不能只看旧 LIVE/status 文件。不要在作业仍活跃时再次提交相同输出任务。

安装完成后不要改动源码、依赖库、模块或 config.sh；全局设置与环境冻结是恢复和跨任务汇总的前提。需要改资源时可以用 `MSI_MEMORY_OVERRIDE=32G MSI_TIME_OVERRIDE=48:00:00 bash msi/submit.sh chains`，无需改研究参数或冻结文件；申请值应在分区限额内。需要升级软件/更换路径时，解压一份新包并建立新的 library/results。

安装失败可诊断日志后使用新的 library 路径重做 setup；不要对有研究结果的安装库原地升级。汇总失败会保留任务原件；完整集合已有时会核对原件 hash 后继续验收。

## 包内文件与当前研究状态

- `IntegMultiReg_0.2.0.tar.gz`：已验证候选源码包，SHA256 `00d557a648f27663996fbd6e9270759d3499b59e693b9388d07268d00c00cec2`。setup 只安装它。
- `dependencies/`：32 个锁定版本的非 base R 依赖源码、按依赖排序的版本表和下载 URL/hash。Matrix 1.7-1 要求 R ≥ 4.4。
- `paper/`、`Ref/`：研究 runner、数据、原始生成器、稿件和参考材料。MSI runner 只新增任务选择/记录；原版三文件及精确 diff 在 `provenance/`。
- `IntegMultiReg/`：含未提交修改的开发源码快照，供后续修改和交接；不是 setup 的安装输入。git HEAD/status/diff 在 `provenance/`。
- `output/`：已完成的原规模八链、Table 1、两类单重复模拟、诊断、最新稿件 PDF、严格 CI 证据以及相关说明。Mac 的绝对路径/hash 保留作历史来源，不是 Linux 的运行依赖。
- `handoff/`、`START-HERE.md`：历史交接；以本文件和最新状态说明为准。

严格 native/Valgrind 验收已通过。原八链已完成，但甲基化 rank Rhat≈1.2007、bulk ESS≈27，混合不足。Table 1 原规模 code2017 运行及验收完成；主模拟/相关模拟只有每配置一次的规模测试，完整 50/30 重复尚未完成。Mac 上的加长第一条链中断，未留下完成 checkpoint，原因没有确认。论文结构/图表已修订，但仍待最终科学结果和文字复核。

不包含 Mac 安装库、编译产物、.git、重复历史大压缩包和原始下载缓存。已有整理数据和生成器足以运行任务；原始数据下载脚本保留。`provenance/inventory.json` 记录打包范围。历史记录提到的部分旧 pilot/截图未全量携带，不影响当前运行。

## MSI 官方参考

- [Slurm 作业提交](https://userdocs.msi.umn.edu/compute/slurm_job_submission.html)
- [共享分区与时限](https://userdocs.msi.umn.edu/compute/shared_partitions.html)
- [计算节点与 login node 的使用](https://userdocs.msi.umn.edu/compute/cluster_info.html)

本地已通过源码安装、真实拆分任务对照和负例测试；具体结果见 `provenance/LOCAL-VALIDATION.md`。完整归档解压/hash 检查另记录在随包提供的 `IntegMultiReg-MSI-20260923.VERIFY.json`。没有访问或提交 MSI，因此真实模块兼容性、调度及 Linux 安装必须由 setup + smoke 验证。
