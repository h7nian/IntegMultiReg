# Anvil 当前工作区

完整上传包已独立解压，并在开发前验证全部 3381 个文件。
当前源码、研究入口和调度入口分别在 `IntegMultiReg/`、`paper/`、`hpc/`。

先读 [当前实施状态](IMPLEMENTATION.md) 和 [Anvil 执行说明](hpc/README.md)。
`SOURCE-VERIFIED.json` 记录恢复来源；Git 基线保留恢复时的源码。
原目录以及冻结历史结果不覆盖，历史 MSI 脚本在 `provenance/msi-before-anvil/`。

依赖构建、候选回归、sanitizer、九项 GitHub CI 和严格 Valgrind 均已通过。
查看 [自动更新的进度与阶段结果](runs/study-20260923/PROGRESS.md)，每 15 分钟刷新。
按用户要求，可用 `--pending-native-job JOB_ID` 让原法审计和完整规模 pilot
与严格检查并行；全量任务和最终汇总仍要求严格检查通过。不能从旧 CI
通过记录或 Slurm 作业结束推断当前候选已通过，也不能手工创建验收标记。

方法参数仍由同一套公开 API 处理：code2017 原法、paper 严格零 ridge 审计，
以及单列的 paper/ridge=0.001 对照。模拟每配置每次 replicate 是一个 task，
当前研究按用户要求将总并发上限提高到 64：剩余试跑最多 1 个、全量先 63 个，
试跑验收通过后全量升至 64。执行策略与控制状态在
`runs/study-20260923/execution/`；科学源码和既定重复数保持冻结。
完整模拟重复数为主模拟 50、相关模拟 30。
