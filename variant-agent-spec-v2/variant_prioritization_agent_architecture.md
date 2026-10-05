# Variant Prioritization Agent — Version 2 Specification Package

版本：2.1.0-draft（review 后修订，见 [REVIEW.md](REVIEW.md) 与 [DECISIONS.md](DECISIONS.md)） · 日期：2026-10-05 · 状态：可交给其他 coding agent 审核的工程规范

这份文件更新原来的单文件架构说明，作为新版阅读入口。核心文件用英文编写，便于直接交给 coding agent；本页说明整体选择和需要你审核的科学问题。它们定义的是待实现系统，不代表代码已经写好、实验结果已经验证，或排名政策已经获得科学认可。

## 1. 系统目标

优先实现：**疾病 + gene + 实验背景 → one genomic allele per row 的候选变异表**，用于 scPRIME 或其他 variant-level functional screen。

疾病→risk gene 是独立的未来轨道；新版保留扩展接口，但不在第一版同时实现。SETD1A、schizophrenia、iPSC-derived neurons、scPRIME 用作首个开发案例，而不是全系统泛化能力的证明。

## 2. 文件分工

| 文件 | 回答的问题 |
|---|---|
| [SPEC.md](SPEC.md) | 系统做什么、科学边界、具体输出、第一版如何排序、什么算完成 |
| [EVIDENCE_SCHEMA.md](EVIDENCE_SCHEMA.md) | Variant、claim、source、dataset、editing、validation 如何存储和关联 |
| [ARCHITECTURE.md](ARCHITECTURE.md) | Workflow、agent loop、tools、state、checkpoint、预算、失败和重放如何运行 |
| [EVAL.md](EVAL.md) | 如何建立 gold、计算指标、避免泄漏、区分工程正确与科学有效 |
| [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) | 12 个可分配开发任务，依赖、验收和第一轮 builder prompt |
| [REVIEW_GUIDE.md](REVIEW_GUIDE.md) | 下一位 reviewer 要挑战哪些假设、怎样报告问题、可直接复制的 review prompt |

科学行为以 SPEC 为准；对象语义以 schema 为准；执行以 architecture 为准；指标分母以 eval 为准。有冲突必须指出并修改，不允许实现者默默选择。

## 3. 新版的关键决定

### LLM 的权限被限定在证据工作里

LLM 可以选择允许的检索、解释杂乱文本、提出 evidence claims、判断 dataset relevance、生成明确标注的机制假设。Host code 决定状态迁移、预算、工具权限、identity、证据验收、排名与导出。

固定工作流里的 semantic worker 不一定需要自主 agent loop。不同角色可以由同一个模型执行；不要求为了架构名称同时运行多个 agent，也不要求使用 LangGraph。

### 精确 allele 是最重要的主键

蛋白简称、rsID、截断 frameshift 描述都可能对应多个 genomic alleles。系统保留这些 mentions，但必须解析到 reference-validated GRCh38 allele 后才能进入正式排名。保留 transcript versions、原始 assembly 和映射依据。

### 证据强度不再只有一个 level

分别表示：是不是这个 exact allele、是不是这个疾病、什么实验、什么细胞背景、什么结果、比较组和质量如何。Gene knockout、same residue、domain、gene burden 都不能自动变成 exact-variant functional evidence。

ClinVar submission 数、提交机构数、文献数和去重后的患者/家系数分开。GEO 必须区分 gene expression measured、gene perturbation、disease dataset、genotype-confirmed exact allele samples。

### 第一版排序：allele nomination（2.1 修订）

目标是为下一轮 screen 提名**新的** allele，并可迁移到其他 gene/disease（DEC-01）。默认 profile `nomination_v0.1`：在 (channel, consequence stratum) 内按 `(K, O, P)` 排序——K 为 gene 机制与变异类别是否吻合（LoF/GoF/dominant-negative/unknown，class-level 推断，不是 allele 证据），O 为患者 cohort 观察，P 为该类别声明的单一预测器原始值。已有 exact-allele 功能实验的 allele 标记为 characterized，默认进入 control channel。并列显式报告，不用 hash 打破。跨类别分配由独立的 panel 步骤完成。旧的 `evidence_priority_v0.1` 只作为 eval baseline。

编辑可行性和生物相关性分开。缺少编辑工具、host genotype 或 editor 参数时输出 not_assessed；不能把 PAM 存在变成预测编辑效率，或称表格为可直接下单的 optimized panel。

### “可重放”与“模型稳定”分开

Replay 固定 accepted evidence 和 annotation，不重新调用网络或模型，要求相同内容与排名。Re-extraction 在相同文档上重新调用模型，测量变化。Refresh 更新来源并保留漂移记录。temperature=0 不是整个系统结果完全相同的保证。

## 4. 推荐推进方式

1. 把整个 ZIP 交给其他 agent，使用 REVIEW_GUIDE 的 reviewer prompt。
2. 让它先交 REVIEW.md，明确 P0/P1 问题和具体修改建议。
3. 审核你关心的科学政策，尤其是排序目标、相关疾病、controls、editing 和 gold。
4. 再分配 ISSUE-001–004，先建立可离线运行、可重放的最小闭环。
5. 验收之后再连接真实数据库和 literature/GEO；没有 gold 时如实输出 not_evaluated。

不需要你先解决所有远期问题，基础工程可以按已写明的 provisional 默认推进。但 coding agent 必须把提出的新科学假设记录下来。

## 5. 最值得你亲自审核的选择

| 决策 | 新版默认 |
|---|---|
| 首要目标 | 新 allele nomination（已确定，DEC-01） |
| 相关疾病 | 展示并标为 related，不自动当成 schizophrenia-specific support |
| Transcript | MANE-oriented 展示，保留其他 isoform effects；identity 仍是 genomic allele |
| LoF/missense/splice 比较 | 保留类别与机制，不用一套未校准 raw prediction score 混排 |
| Controls | 先采用明确提供的候选与依据，不把 synonymous 自动视为 benign |
| 编辑可行性 | 有实际工具/设计/实验才给相应结果，否则 not_assessed |
| Gold | 需要你提供或审核 exact alleles、来源和 label purpose；不根据记忆编造 |

## 6. 相比旧版的变化

保留：agent/deterministic 分界、evidence store、来源追踪、实验可行性、长期 risk-gene 轨道。

新增或收紧：可验收 requirement/invariant IDs、reference-aware identity、claim/link 分离、orthogonal evidence dimensions、candidate expansion closure、拒绝证据审计、sample-level GEO genotype、source capability gates、严格 replay、evaluation denominators、leakage controls、分 issue build 和 reviewer 任务。

旧版泛化的数字总分、信息增益公式和 assay 等级被替换为明确标注的待验证政策；它们可以作为未来候选方案，不能在第一版被实现成已经验证的科学事实。

## 7. 完成状态

已完成：新版规范与审核/开发任务定义。

尚未完成：外部 agent 审核、代码实现、source live integration、gold adjudication、真实 benchmark、editing engine 选择与实验面板优化。下一步 review 可以直接使用这套文件，不需要依赖聊天记录。
