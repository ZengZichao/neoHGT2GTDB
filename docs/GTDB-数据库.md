> 语言:中文(本文件,与英文版内容对照)| English (primary): [docs/en/GTDB-database.md](en/GTDB-database.md)

# GTDB 数据库（实验性）

neoHGT2GTDB 支持以 [GTDB](https://gtdb.ecogenomic.org/) 分类体系构建参考库 —— 这是本适配对上游
[neoHGT2](https://github.com/cgneo/neoHGT)（命令行入口同为 `neoHGT`）的深度改进所在。GTDB 提供标准化、
层级完整的原核分类（每个物种簇恰好一个代表基因组），可显著改善 self/close/distal 分组的
统计性质。GTDB **不含真核与病毒**，需要时使用 `--source hybrid`（该模式需要另行获取 NCBI
taxdump，见"已知限制"一节）。

本页的所有数字均为本机实测（GTDB R09-RS220 本地数据集，2026-11-01/02）。
英文对照:完整版见 [`docs/en/GTDB-database.md`](en/GTDB-database.md),README 中另有 "GTDB database (experimental)" 英文精简节。

## 工作区约定（必须先建立）

所有产物只能写入三个目录，且 `--gtdb-dir` 指向的数据集**全程只读**（工具会在启动时自检写
权限，并拒绝任何指向它的输出路径）：

| 目录 | 用途 |
| --- | --- |
| `neoHGT2GTDB-项目代码/` | 源码、测试、文档 |
| `neoHGT2GTDB-临时文件/` | `--tmpdir`、DIAMOND 中间态、可丢弃 |
| `neoHGT2GTDB-缓存测试/` | `db/`、`search/`、`analyze/`、`fixtures/`、`logs/`、`benchmarks/`、`reports/` |

先加载环境（会重定向 matplotlib/XDG/TMPDIR 到项目内，避免写到工作区之外）：

```sh
. tools/env.sh
```

不过 GTDB/hybrid 建库**不再依赖这个 wrapper**：上游在没给 `--tmpdir` 时把暂存目录开在
`mkdtemp()`（即环境 TMPDIR 所指处），所以以前"中间产物在工作区内"其实是 env.sh 的副作用。
现在 `--source gtdb|hybrid` 分支会把它迁到 `workspace.tmp_dir` 下（`临时文件/gtdb_XXXX`），
并顺手 `rmdir` 掉上游刚建好的那个空目录（`rmdir` 只能删空目录，删不到别的东西）；构建结束时
照旧整目录清理。实测：把 `TMPDIR` 指向别处、不给 `--tmpdir` 跑一次 1 基因组构建，日志打印
`Scratch directory inside the workspace: …/neoHGT2GTDB-临时文件/gtdb_7sh52ard`，被当作诱饵的
那个目录**始终是空的**，构建结束后也不留残留（`缓存测试/db/nomprefail`）。NCBI 路线的默认值
保持上游原样（RT-2：不改既有选项默认行为），所以走 NCBI 建库仍应先 source env.sh。

三条约定可以用 `tools/verify_workspace.py` 机械核对（IT-7 / DoD-9）：

```sh
python tools/verify_workspace.py --out "$NEOHGT_CACHE/reports/IT7_workspace.json" \
    --since "2026-11-01 17:25:00" --probe-report \
    "$NEOHGT_CACHE/reports/IT7_tmp_disposable_probe.json"
```

它检查四件事：(1) 三个项目目录之外没有**本次项目产生**的文件（早于起始时刻的既有文件
单列 `predates_project`，不算违规）；(2) 项目目录内没有指向工作区之外的符号链接；
(3) 缓存目录保持上表的子目录结构；(4) 数据集 mtime 审计（复用 `verify_read_only.py`）。
"**临时文件可安全清空**"这一条不是靠读清单猜的：清单里出现的暂存路径只作为**清单**列出，
真正的判定靠行为探针——用一个**全新的空 `--tmpdir`** 续跑构建，再由
`tools/verify_resume.py` 与首次构建逐文件比对（PASS 即证明旧的暂存内容无人依赖）。

## 只读预检（建议每次构建前执行）

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/logs/layout_probe.json"
# 等价入口（不经 CLI 包装）：
python -m neoHGT.gtdb probe --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

两种入口等价：`neoHGT gtdb <子命令>` 只是把参数原样转交给 `neoHGT.gtdb.main()`
（含退出码；命令行入口 `neoHGT/gtdb.py::cli_main` 会把"输入不合格"这类
`ValueError`/`SelfCheckError` 打成一行 `ERROR: …` 并返回退出码 2，而不是甩一段
回溯——`main()` 本身照旧抛出，库调用方不受影响），子命令为
`probe / taxdump / map-report / gtdbtk / reconcile / tree-dist`；更早的扁平写法
`neoHGT gtdb --action probe|taxdump|map-report` 仍可用。

输出布局标签、release、序列 ID 约定与锚点覆盖率，并执行 **19 项**自洽核对
（faa 文件数 vs metadata 代表集数、taxonomy 行数 vs metadata 行数、谱系深度恒 7、
rank 前缀词表、六档 `*_Sample_<rank>.IDlist` 行数 vs 唯一类群数、域集合、
accession 对齐、代表基因组是否有本地序列）。**任一不符即中止构建**（可显式
`--gtdb-no-selfcheck` 跳过，但结果不再保证可复现）。

本地 r220 实测：布局 `gtdb-r220-renamed`，`_SerialNumber_` 锚点覆盖率 1.00，
最长 defline 270 字符，19/19 核对通过，耗时 28 s。

## 不符清单落盘（DoD-2）

`probe` 只给出**计数**是否相符；要把"哪几个 accession 不符"落到磁盘上，用
`reconcile`：

```sh
python -m neoHGT.gtdb reconcile --dir "$GTDB_DIR" \
    --out "$NEOHGT_CACHE/reconcile" --report "$NEOHGT_CACHE/reports/DoD2_reconcile.json"
```

它在 metadata 全表、GTDB taxonomy 与本地 `.faa` 三方之间做集合差，写出五个
TSV（每个都有一行 `#` 注释头说明来源/ release / 生成时间与一行 `accession`
表头，因此**即使为空也留痕**），不一致时退出码为 1：

| 文件 | 含义 |
| --- | --- |
| `missing_faa.tsv` | metadata 标记为代表、但本地没有序列文件的 accession |
| `extra_faa.tsv` | 本地有 `.faa`、但 GTDB metadata 完全不认识的 accession |
| `nonrep_faa.tsv` | 本地落盘的**非代表**（成员）基因组文件，默认布局下应为空 |
| `metadata_not_in_taxonomy.tsv` | metadata 有、`*_taxonomy.tsv` 没有 |
| `taxonomy_not_in_metadata.tsv` | taxonomy 有、metadata 没有 |

本地 r220 实测（`reports/DoD2_reconcile.json`）：metadata 596,859 基因组
（其中代表 113,104）/ 本地 FASTA 113,104 / taxonomy 596,859 行，**五张清单
全为 0**，`consistent: true`；耗时 24.3 s，峰值内存 1.7 GB。

## 单独构建 taxdump

```sh
python -m neoHGT.gtdb taxdump --dir "$GTDB_DIR" --out "$NEOHGT_CACHE/db/r220_taxdump"
```

实测：**146,380 节点**（GTDB 146,379 + root），rank 计数 = domain 2 / phylum 217 /
class 607 / order 2,030 / family 5,460 / genus 24,959 / species 113,104；
转换 42.9 s，峰值内存 1.5 GB；产物 192 MB。

**只有谱系文件的数据集也能出 taxdump**。分类学树本身是谱系串的函数，不需要
metadata 表；官方 release 解包后把 `bac120/ar53_taxonomy_rNNN.tsv(.gz)` 直接放在数据根目录，
原实现的发现逻辑只往 `tree_bac120/`、`tree_ar53/` 里看，于是"官方布局"反而不被识别。现在：

```sh
# 真实数据实测：本地 r214 副本（无 metadata 表、谱系文件在根目录）
neoHGT gtdb taxdump --dir /path/to/GTDB_r214/parent_directory \
    --gtdb-id-manifest <r220 库>/gtdb/taxid.map.gz --out <新目录>
```

打印 `NOTICE: no metadata table ...` 与 `lineage_only: true`，只跑"没有 metadata 也能核"的那一项
（谱系深度恒 7、同一 accession 不得有两个谱系），并明确列出被跳过的 6 项 §3.6 核对；实测
113,318 节点 / 6.3 秒 / `official-plain`，同时写出 `gtdb/{taxid.map.gz,delnodes.txt,reserved_ids.txt}`
（`缓存测试/gtdb_r214_taxdump_v2/`；第一次跑的 `缓存测试/gtdb_r214_taxdump/` 里那份"被跳过"清单**名字是编的**，
节点数与映射仍有效——见本文件"自我发现"第 10 条）。
建库（`neoHGT database --source gtdb`）仍然**要求** metadata 表，报错会点名 `--gtdb-metadata`。

谱系串不合格时在**任何写入之前**被拒绝——本地 r226 副本正是这种输入（第三方把等级用 `_` 粘连、
并截断到属级，accession 列还残留 `.faa` 后缀）：

```
[FAIL] taxonomy_files_readable: .../ar53_taxonomy_r226.tsv:1: lineage of RS_GCF_945873965.1.faa
       does not have 7 ranks (provenance check failed).
```
退出码 1，输出目录不会被创建。核对与双向（旧映射接新谱系、旧映射接更旧谱系）的实测数字见
`缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`。

## 建库

```sh
# 极小档实测（r220 古菌每科 1 个代表，564 基因组 / 905,829 蛋白，墙钟 2 分 37 秒）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" --rank-sample family \
  --domain archaea --compile diamond -o "$NEOHGT_CACHE/db/r220_ar_family"

# 谱系子集（示例：肠杆菌科代表集）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --seq-source local \
  --compile diamond --threads 4 -o "$NEOHGT_CACHE/db/r220_entero" \
  --tmpdir "$NEOHGT_TMP/build/$$"

# 使用外部基因组清单（可直接吃本地现成 *_Sample_<rank>.IDlist）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --genome-list "$GTDB_DIR/metadata_bac120/bac120_metadata_r220.tsv_Sample_Genus.IDlist" \
  --compile diamond -o "$NEOHGT_CACHE/db/r220_bac_genus"

# 全量代表集（113,104 基因组 / 333,384,550 蛋白，零下载）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" --reps --seq-source local \
  --min-completeness 50 --max-contamination 10 \
  --compile diamond --keep-faa compress \
  -o "$NEOHGT_CACHE/db/r220_reps" --tmpdir "$NEOHGT_TMP/build/$$"
```

`--max-genomes N` 会截取选集前 N 个基因组，用于受控冒烟测试。
质量类筛选走 `--min-completeness / --max-contamination / --mimag / --type-material gtdb|ncbi`
（上游的 `--typemater` 在 GTDB 路线就是 `--type-material ncbi` 的别名；`--complete` 会被**拒绝**而不是
静默失效——各列的真实计数见下面"质量类选项在 GTDB 下的语义"一节）。
`--gtdb-dir / --seq-source / --local-layout / --gtdb-taxid-base` 未显式给出时取默认——`--seq-source` 默认 `local`、
`config.yml` 的 `database.gtdb`（`dir` 即本地只读数据集根），命令行始终优先。

抽样三件套 `-s/--sample`、`-r/--rank`、`--above` 在 GTDB 路线的语义与 NCBI 路线一致
（本次更新补齐：此前 `--above` 没有进入 provider 的 opts，GTDB 下被静默忽略）：

- 先把每个 `<rank>` 类群取至多 n 个基因组；
- 谱系串**没有到达该秩**的基因组（例如只鉴定到属 `g__G2` 而没有 `s__`）在这一轮
  不参与，程序会打印"多少条未鉴定到该秩"并提示 `--above`；
- 给了 `--above` 之后，从采样秩的上一级开始逐级（至门级为止，与上游一致不把域当作抽样档）
  把每个粗类群补足到 n 个——这正是上面那类"浅谱系"基因组进入库的通道。
- `--rank-sample <rank>` 是另一条低成本档：每类群只取质量分最高的 1 个，不看 `--above`。

一个诚实的实测结论：r220 的**代表集谱系一律 7 段到种**（`genomes.tsv` 里
`gtdb_taxonomy` 的深度直方图是 `{7: N}`），而"细类群每档取 n 个"必然让每个粗类群
也至少满足 n（除非候选已经取尽），所以在全深度数据上 `--above` 是**可证明的 no-op**：
同一家族同一组参数实测 235 vs 235 个基因组。它起作用的场景是谱系更浅的输入
（GTDB-Tk 结果、成员基因组扩量），这条路径由单测锁定
（`test_above_tops_up_coarser_groups`：一个只鉴定到属的类群只有 `--above` 才进得来）。

`search` 端已把 `--keep-placeholder-names` 接上（默认 `yes`，来自
`config.yml` 的 `taxonomy.profiles.gtdb.keep_placeholder`）：GTDB 里有 1,658 个类群的
名字是编号（`2-12-FULL-45-22`、`g__DUMJ01`、`s__MK2026 sp010841923`），它们**是**官方
标签；如果用户显式打开 `--tax-capital` / `--tax-latin`，`util.is_placeholder_gtdb()`
会豁免这类标签，否则上万个命中会静默消失。`--keep-placeholder-names no` 回到严格语义。
豁免只在 GTDB 数据库上生效，且只在那一次过滤调用期间替换判定函数，调用结束即还原。

> 这段语义有两处修正（原实现在这两种情形下都是坏的）：
> - **`--tax-latin yes` 在 GTDB 下过去必然丢掉全部命中**。上游 `taxfilt_wf()` 先解析出种级
>   **taxId**，再把这个 taxId 当作"名字"交给 `is_latin()` —— 数字串永远不是二名法拉丁名，于是
>   每一条命中都被判为"非拉丁"。GTDB 路线现在先把该 taxId 换成它自己的名字再判定
>   （`neoHGT/search.py::taxfilt_wf_gtdb`）；NCBI 路线一字未动（RT-1 逐条比对通过），
>   这一上游缺陷记在 内部QC台账《计划任务对照表》 的 3.5 行。
> - **严格语义不等于拒绝带前缀的官方名**。`--keep-placeholder-names no` 现在只是取消编号标签豁免，
>   判定仍走前缀容忍版（`s__Halobacterium` 依然算大写），因为 DoD-5 要求的是"按去前缀后的裸名判断"。
>   用例：`neoHGT/tests/test_gtdb_pipeline.py::TestSearchChoiceValues` 与
>   `TestPlaceholderExemption`。
>
> 同一轮还修了 `--require-staxids`：它被解析、被 `arg2bool` 转换、被存进对象，**但没有任何代码读它**，
> 所以 `no` 完全无效；而且因为 `arg2bool(None)` 是 `False`，"用户没提这个选项"会被当成"用户不要"。
> 现在 GTDB/hybrid 路线默认 `--require-staxids yes`：命中表缺 staxids 列立刻报错（RK5 的保护）；
> 只有显式 `--require-staxids no` 才降级为一条 WARNING 并继续（此时分类只能靠 subject Id）；
> NCBI 路线恢复上游语义（缺列沿用上游行为），该默认不再作用于 NCBI 运行。

`--compile none` 也修了一处：写完 `db.faa`/taxdump/侧表之后 `write_gtdb_meta()`
把还没定位过的 `self.diamond=None` 交给 `shutil.which()`，抛 `TypeError` 让整个构建在
最后一步失败——看起来"库已经建好了"却没有 `meta.json`。现在先判空，`meta.json` 的
`diamond_version` 记为 `null`，并由 `TestMetaWithoutDiamond` 锁定（含"DIAMOND 存在时
仍然如实记录版本"的另一半）。

方案 §6 的编号任务与实现的逐条对应关系写在
内部QC台账《计划任务对照表》(工作区侧,未随公开版发布)，并由可追溯门禁核对（编号齐备、
状态取固定词表、引用的路径确实存在、"部分实现/未做"必须写清缺什么）。

方案附录 F 列的 14 条"与本适配同源的既有缺陷"登记在
内部QC台账《既有缺陷登记》(工作区侧,未随公开版发布)：每条要么是**已修复**（附证据），要么是
**登记未改**并写清为什么不改（多数是"修它必须删除上游已有行"，与 RT-1 的纯追加约束冲突）
以及如何复现。同一个门禁也核对这张表，且行号引用同样会被验证（`search.py:1166` 这种写法
如果指向不存在的行会被报出来）。

方案的 §9 风险登记册（19 条）核对在
内部QC台账《风险应对登记》(工作区侧,未随公开版发布)：每条写"计划里的缓解最后落在哪里"，状态取
`已实现缓解 / 部分实现 / 登记接受`；今天的成绩是 16 已实现缓解、3 部分实现
（RK7 的已知 HGT 灵敏度验证、RK9 的远端断点续传执行、RK19 的"NFS→NFS 输出"备选，
分别卡在 NCBI 参考库授权、联网授权与"不得写工作区之外"的安全规则）。

第八、第九张：方案 §2.3 对 `search.py` 逐行写下的 9 条 NCBI 假设核对在
内部QC台账《搜索层假设核对》(工作区侧,未随公开版发布)（5 条已按 GTDB 改造、3 条确认无需改动、
1 条**与方案当时的结论相反**——方案认为 BLAST 路线"可行"，实测 `makeblastdb` 无法为
151–270 字符的 GTDB 序列 ID 建索引，于是 `--compile blast/both` 改为默认拒绝、`--force`
才继续，见 ADR-8）；方案 §0.5 的 19 条审查修订核对在
内部QC台账《方案修订记录核对》(工作区侧,未随公开版发布)（18 已采纳、1 部分采纳：REV-12 的夹具来源与
§10.1 覆盖点⑥⑦未照搬，逐条写在交付物文档的核对表里）。这两族同样是"带编号、可判真伪"的
方案内容，此前没有任何登记处。

第十张是 §13 里程碑核对（内部QC台账《里程碑核对》）。
在这些表之外，门禁本身现在也进了 CI：`docs/方案编号清单.json` 是把方案里全部编号清单
（276 项：28 任务 / 14 缺陷 / 12 交付物 / 19 风险 / 10 DoD / 11 决策 / 36 验证矩阵 /
9 搜索层假设 / 19 审查修订 / 37 测试编号 / 27 §7 选项 / 27 附录 C/D 选项 / 27 附录 B 符号）冻结下来的摘要，
所以**方案文档不在仓库里也能核对**。CI 的 "Register gate" 步骤用这份摘要跑十张表；
只有在本地拿着真方案时才可能重算摘要（`--plan <方案> --emit-ids docs/方案编号清单.json`），
而用例 `test_frozen_plan_digest_equals_the_plan_it_came_from` 会在本地断言摘要与方案逐字段相同，
摘要一旦过期就红。另外这个工具自己也拒绝把工作区之外的路径当输出目标
（`guarded_out()`，规则 11：`--out /tmp/x` 不是重定向，不受 `tools/run.sh` 拦截，所以得自己拦）。

方案的 §10.2 / §10.3 / §10.4 / §10.5 验证矩阵（36 条编号：20 个 UT、9 个 IT、4 个科学有效性项、3 个 RT）核对在
内部QC台账《测试用例对照》(工作区侧,未随公开版发布)：每条写明由哪个用例（带 `文件:行号`）或哪份证据产物
覆盖，状态取 `已覆盖 / 部分覆盖 / 未覆盖`。当前 33 已覆盖、0 部分覆盖、3 未覆盖（IT-3 的补数路线、10.4-1 与 10.4-4 的 NCBI 参考库正面对照，都卡在联网未获授权）。§10.4-3 的参数敏感性已有实测曲线：close_size、input_cov、bandwidth、noise 四个旋钮都有实测曲线，见 `缓存测试/reports/S10_4_3_sensitivity.md`。这张表源自一次清点：当时发现
§10.1 要求的派生夹具并不存在，而 `search`/`analyze` 的 GTDB 分支从未有过**离线**端到端用例——
把验证矩阵摊开成可机器核对的表，就是为了不让"计划里写了、代码里其实没有"再静默发生。

其中 RK11 本次更新补齐：GTDB 会换掉高阶分类单元的名字（`p__Proteobacteria` →
`p__Pseudomonadota`），而官方 `synonyms.*.tsv` 不在本地副本里，所以
`analyze` 新增 `--taxon-synonyms <alias<TAB>current 的 TSV>`；`--self-tax`、
`--close-tax`、`--tax-include`、`--tax-exclude`、`--input-tax` 五种引用都会先过
别名表再解析，用到旧名时打印一条"这是弃用名，已改用 X"的提示。

方案 §1.4/§10.7 的 10 条验收标准（DoD-0…DoD-9）对照在
内部QC台账《验收标准对照》(工作区侧,未随公开版发布)：状态取 `达标 / 部分达标 / 被外部条件阻塞`，
每条给证据路径。现状 9 达标 + 1 部分达标（DoD-2：代表集**基因组数**由构建期强制项对账，
**蛋白条目数**是 metadata 求和 + 分档实测外推，因为全量抽取按使用方决定未跑）。
同一轮把 §3.6 的三行全库文件数清点（`.faa`/`.fna`/`.tsv`）做进了只读审计：
`tools/verify_read_only.py` 顺带输出 `files_by_extension`，实测
`faa=146686 fna=113360 tsv=112972`，与 §3.6 表内数字逐项相等
（`缓存测试/logs/read_only_check_r12.txt`）。这三行取决于副本的具体内容，所以按
"每次审计留痕"而不是"构建期等值断言"处理，`--no-inventory` 可关掉这次清点。

方案 §8 的 11 条设计决策（ADR-1…ADR-10 含 4b）落在
内部QC台账《决策落点对照》(工作区侧,未随公开版发布)：状态取
`已落实 / 部分落实 / 偏离并说明 / 未落实`，逐条给代码落点。现状 10 已落实 + 1 部分落实
（ADR-4 的"v1 只承诺代表集"在代码里成立，但其可跑通性证据只到 3.25M 蛋白档，全量档未跑）。

方案的 §12 交付物（D1–D12）与**附录 B/C/D 的接口形态**核对在
内部QC台账《交付物清单》(工作区侧,未随公开版发布)：12 项交付物逐条给状态与产物位置（D11 全量档基准
与 D12 的 NCBI 正面对照如实标为部分交付），附录 B 承诺的 27 个函数/类/常量必须能在代码里
找到（允许落在同包另一个模块），找不到就必须在该文件的"附录 B 符号偏差"表里登记并说明；
附录 C/D 示例命令里的 27 个 `--flag` 必须要么是本工具 CLI 的真实选项、要么是第三方工具的
选项、要么记在"附录 C/D 命令形态偏差"表里。

> 上表第二条（肠杆菌科）为实测命令：`f__Enterobacteriaceae` 在 r220 恰有 **793** 个代表
> 基因组（153 属），`genomes.tsv` 亦为 793 行，与 metadata 逐一相符；未叠加质量阈值
> （`--min-completeness/--max-contamination` 需显式给出，见上）。

> **质量类选项在 GTDB 下的语义**：本地 r220 metadata 全表 **596,859** 行
> （bacteria 584,382 / archaea 12,477），其中代表集 **113,104** 行。下面的计数全部由**生产选择代码**
> `neoHGT/gtdb.py::GtdbProvider.select()` 在这张全表上实测得到（`缓存测试/reports/R29_quality_class_counts.json`），
> 不是抽样估计——首次记录时只读了表的前 40 万行，得到"带 designation 的只有 14,745 条 / 3.7%"，
> 那是**取样偏差**（type material 集中在文件后段），已用全表数字替换。
>
> - `--type-material gtdb` → `gtdb_type_species_of_genus == 't'`：全表 **4,545** 行，代表集 **3,347** 行。
> - `--type-material ncbi`（上游旗标 `--typemater` 在 GTDB 模式就是它的别名）→ `ncbi_type_material_designation`
>   且必须排除 `na`/空串：全表 **22,256** 行（`na` 574,603），代表集 **16,371** 行＝14.5%，分 5 种
>   designation：type material 16,313 / synonym type material 28 / reftype 19 / neotype 9 / pathotype 2。
> - `--mimag high|medium|low|any` → `mimag_*_quality == 't'`：代表集 27,292 / 85,752 / 8，三档互斥，
>   `any`（§7 的第四个 choice，含义是"必须带有某一档 MIMAG 分级"）在代表集留 **113,052 / 113,104**、
>   全表留 596,497 / 596,859，即"没有任何分级"的基因组只有 52 / 362 个。列缺失时报错而不是静默返回全表：
>   GTDB-Tk 输入表和不发布 MIMAG 类的 release 就没有这些列，此前 `--mimag` 在那种表上是
>   **静默 no-op**。
> - `--complete` 在 GTDB 下**直接报错**而不是静默失效：上游看的是 `ncbi_assembly_level`，r220 代表集上它
>   只留下 8,593 个（Chromosome 1,070 + Complete Genome 7,523），被丢掉的 104,511 个中位 completeness 仍有
>   94–96%——assembly level 描述的是"序列怎么提交的"，不是"基因组完整不完整"（`ncbi_refseq_category`
>   同样不能用来判质量：代表集 87.4% 是 `na`）。方案 §7 R5 明确以 `--min-completeness/--max-contamination`
>   替代，实测在代表集上 completeness≥50 留 82,662、contamination≤10 留 82,464、两条同时 82,161。
> - 构建报告另起一行输出这些类的**实际计数**（`quality classes in the dataset: …`），因为方案 §7 要求
>   报告给出真实数字而不是让使用方猜；计数在**未过滤**的表上统计，所以它描述的是数据集而不是本次筛选。

### 库目录结构

```
<db>/db.faa  genomes.tsv  lineages.txt  taxon.map.gz  genome.map.gz  meta.json
<db>/taxdump/{nodes.dmp,names.dmp}     # 全量 GTDB 树（rank 用 NCBI 词表）
<db>/diamond/db.dmnd
<db>/gtdb/taxid.map.gz                 # 跨版本 taxId 稳定性凭据（只追加不重排）
<db>/gtdb/{labels.tsv,lineages.tsv,genome2taxid.tsv.gz,stem2taxid.tsv.gz,
          prefix2genome.tsv.gz,ncbi_map.tsv.gz,build_report.md,
          input_fingerprints.tsv,extraction_state.json,fetch_plan.tsv,
          members.tsv.gz,delnodes.txt,reserved_ids.txt}
```

`fetch_plan.tsv` 只在 `--seq-source ncbi|gtdb-tar|ena` 时出现；`members.tsv.gz` 只在
`--expand-members yes` 时出现；`delnodes.txt` 与 `reserved_ids.txt` 只在"复用上版
manifest 且确有节点消失"时出现。

`meta.json` 记录 GTDB release、`local_layout`、`id_convention`、输入 provenance（含
本地改造脚本清单与 sha256 指纹）、DIAMOND 版本、taxonmap 摘要、19 项自洽核对结果、
统计量（基因组/蛋白/节点/db.faa 体积/最长序列 ID）。

**派生文件的 provenance 头**：`gtdb/` 下每张表与 `lineages.txt`、
`gtdb/input_fingerprints.tsv`、`gtdb/fetch_plan.tsv`、`gtdb/members.tsv.gz`、
`gtdb/delnodes.txt`、`gtdb/reserved_ids.txt`、以及 `reconcile` 的不符清单，第一行都是

```
# neoHGT version=2.0b3; gtdb_release=220; layout=gtdb-r220-renamed; generated=<本地时间>[; 附加字段]
```

因此任何一份产物都能反查"由哪个版本的代码、哪个 release、在什么时候生成"。
两点例外是**故意的**：`taxon.map.gz` 与 `genome.map.gz` 会被 DIAMOND、BLAST 和
`read_prot2taxid` 原样消费，`--input-tax` 文件（`gtdb gtdbtk` 的输出）会被
`dict_from_param` 原样解析，所以它们不能带注释行。读取 `gtdb/` 内部表的函数
（`load_labels`、`load_manifest`、`load_genome2taxid`、`load_stem2taxid`、
`load_prefix2genome`、`load_ncbi_map`）一律跳过 `#` 行；`load_genome2taxid`
现在会显式找到列头行并在缺列时报明确错误。

### 中断续跑（IT-6）

抽取阶段每次成功后会把计数器写入 `<db>/gtdb/extraction_state.json`。再次以**相同命令**
运行到同一 `-o` 目录时：

- `db.faa` 存在且状态文件存在**且选择指纹一致** → **跳过抽取**，并从 `taxon.map.gz` /
  `genome.map.gz`
  重建 DIAMOND taxonmap（两列同值）、BLAST seq2taxid 与 genome.map；
  `genomes.tsv` / `meta.json` 的蛋白数、去重数、无 taxId 数与一次成功构建完全一致。
- 状态文件记录本次基因组选择的 sha256 指纹；**选择变了就不会复用**旧的 `db.faa`
  （否则会写出与 `genomes.tsv` 不一致的库），程序自动改为重新抽取。
- `db.faa` 存在但 `taxon.map.gz` 缺失 → 直接报错并提示 `--overwrite` 重建，
  绝不发布与 `db.faa` 不一致的元数据。
- 需要重抽时用 `--overwrite`（会删除旧 `db.faa` 并重新流式抽取）。

比对"一次成功构建"与"续跑构建"用 `tools/verify_resume.py`：逐文件 sha256 相同记
`identical`；`taxon.map.gz`/`genome.map.gz` 这类 gz 产物去注释后逐行相同记
`text_equal`；只有 provenance 时间戳行不同的文本表记 `timestamp_only`；
`meta.json` 只允许 `build_started/build_finished/extraction_resumed` 与命令行回显的
`tmpdir`（一个位置，不是构建内容）四项不同。比较是**流式**的——`gtdb/lineages.tsv`
在 6 基因组档就有 96 MB，全量档更大，不能整文件读进内存。任何其它差异都是
`mismatched`，判定 FAIL。

## meta.json 字段规范与 GTDB 版本升级 SOP

### meta.json 字段

| 字段 | 内容 |
| --- | --- |
| `tool` / `source` / `seq_source` | 版本串、数据源（`gtdb`/`hybrid`）、序列通道（`local`/`ncbi`/`gtdb-tar`/`ena`） |
| `gtdb_release` / `local_layout` | 释放版本号与布局判定结果（`gtdb-r220-renamed` / `official-plain` / `official-gz`） |
| `id_convention` / `taxonomy.id_base` / `taxonomy.rank_counts` / `taxonomy.n_superseded` | 合成 taxId 规则、起始 base、各 rank 节点数、被取代节点数 |
| `extraction_resumed` | 本次运行是否复用了已有 `db.faa`（IT-6 溯源） |
| `input_provenance` | 只读输入根目录、`source_writable_by_user`（`os.access(W_OK)` 结果）与 `read_only_policy`（说明强制点是 `guard_write` + 构建后 mtime 审计，而不是文件权限）、taxonomy/metadata 文件绝对路径、序列目录、本地改造第三方脚本清单、缺失的官方辅助文件、`_SerialNumber_` 锚点覆盖率、`input_sha256`（文件 → sha256） |
| `counts` | `genomes{,_planned,_missing}`、`proteins`、`proteins_planned`（本次所选基因组按 metadata 应有蛋白数）、`proteins_extracted_vs_planned`（`ok` / `N short` / `not extracted yet`，由 `neoHGT/gtdb.py::proteins_reconciliation` 包装 `neoHGT/gtdb.py::reconcile_protein_counts` 得出，两个字段同源于一个 dict——早期版本按位置解包，把抽出来的条数写进了"应有"字段；基准是**本次所选基因组**（由 `neoHGT/database.py::write_gtdb_meta` 以 `selection=df` 显式传入），不是整个 release 的表）、`proteins_duplicated`、`proteins_without_taxid`、`taxa`、`db_faa_bytes`、`max_seqid_len`、`by_domain_selected`、`representatives_in_dataset`（整个数据集口径，用于 DoD-2 三方计数） |
| `diamond_version` / `taxid_map_sha256` | 编译所用 DIAMOND 版本、taxonmap 摘要 |
| `selection` | 本次选择方式（`--reps`/`--genome-list`/`--rank-sample`/`--gtdb-taxa`/质量阈值）与实际命中数 |
| `selfcheck` / `selfcheck_items` | 19 项核对的总体状态与逐项 `expected/observed/status` |
| `config_snapshot` / `build_started` / `build_finished` | 生效配置快照与时间戳 |

除 `build_started/build_finished/extraction_resumed` 外，两次相同构建的 `meta.json`
必须逐字段相同（`tools/verify_resume.py` 强制检查）。

### 升级到新 release（例：r220 → r232）

1. **先探测**：`neoHGT gtdb probe --dir <新数据集> --gtdb-release 232 --out ...`。
   19 项核对不硬编码 r220 数字（IDlist 计数由当前数据集现算），因此新释放应全绿；
   若 `idlist_vs_taxa[*]` 或 `rep_genomes_have_local_fasta` 失败，说明该副本不完整，**停手**。
   随后跑一次 `neoHGT gtdb reconcile --dir <新数据集> --out ... --report ...`
   （见"不符清单落盘"一节），把差到的 accession 逐条留下磁盘证据，再决定是否建库。
2. **复用旧 ID 清单**，保证 taxId 跨版本稳定（只追加、不重排）：
   ```sh
   neoHGT gtdb taxdump --dir <新数据集> --gtdb-release 232 \
     --gtdb-id-manifest <旧库>/gtdb/taxid.map.gz --out <新库目录>
   ```
   该命令只做分类学侧产物（`<新库>/taxdump/`、`<新库>/gtdb/`（含新的
   `taxid.map.gz`）、`<新库>/meta.json`），不触碰序列；r220 实测耗时 42.9 s、1.5 GB 内存。
   `--gtdb-taxid-base` 必须与旧库一致（默认 `100000000`），否则新旧 ID 空间会重叠。

   **被删除的类群其 taxId 永久保留**：上一版 `gtdb/taxid.map.gz` 旁边的
   `gtdb/reserved_ids.txt` 会被自动读回，新 release 里新增的类群**不会**顶替历史上
   用过的 ID（否则旧的 `scores.tsv`/`hgts/*.txt` 会在不改一个字的情况下指向另一个生物，
   即方案 RK6）。同 release 重建（无增删）时该文件不会被创建，映射逐条不变——
   r220 实测：`146,379` 节点、两次构建的 `lineage→taxId` 全等。

   **这条 SOP 已在真实的跨 release 场景上跑过**（本地 r214/r220 真实数据，无网络）：
   拿 r220 的 `taxid.map.gz` 去接 r214 的谱系，89,172 个共有节点 ID **零改动**，24,145 个
   r214 独有节点与旧表已知 ID **零冲突**、也没有顶替 57,207 个消失节点中任何一个的 ID，
   这 57,207 个节点全部落进 `delnodes.txt`；反向（r214 映射 + r220 谱系）同样零改动 / 零冲突 /
   零复用，`delnodes.txt` 为 24,145 条。判据要留意：分配按 rank 分块
   （`slot = id_base + rank_ordinal * 1e6`），所以新节点落在自己那一段里是**正确**的，
   "新 ID 必须大于旧最大值"是错误判据（第一版核对就误用过它，见
   `缓存测试/reports/ACCEPTANCE.md` §37.2）。
   数字来自 `缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`，
   离线用例 `neoHGT/tests/test_gtdb_core.py::test_release_swap_never_reissues_or_reuses_a_taxid`
   把同样的不变量固定在 CI 里。
3. **正式建库**同样带 `--gtdb-id-manifest <旧库>/gtdb/taxid.map.gz`；建完后核对
   `meta.json` 的 `taxonomy.rank_counts` 与 `n_superseded`，并抽查若干旧节点 taxId 未变。
4. **重建搜索结果口径**：谱系标签会变（GTDB 改名很常见），历史 `hgts/*.txt` 中的
   `g__` 标签需按新 `gtdb/labels.tsv` 重新解释；`gtdb/ncbi_map.tsv.gz` 也要重新生成。
5. **重新标定分组参数**：库规模变化会移动 `closesize` 的可用区间（见"分析"一节实测：
   564 基因组档只标定到 close_size=100——10/50/100 → 371/371/132 条，200/400 在该档
   未单独实测，按同节结论需接近全量代表集的库密度才有意义）。换库后应重跑
   `tools/run_hgt_demo.py --close-sizes 10,50,100,200,400` 再定默认值。
6. **保留旧库目录**（不删除），新结果目录记录 `db 路径 + meta.json sha256` 以便复现。

## 搜索

DIAMOND 库已内嵌 taxId，因此 `staxids` 直接可用；推荐沿用"手工搜索 + 预计算结果"路径：

```sh
diamond blastp --query q.faa --db "$DB/diamond/db" --threads 8 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive > hits.tsv
neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/search/sample" --tmpdir "$NEOHGT_TMP/search/$$"
```

GTDB 模式下程序会**强制 `--fetch-enable no`**（合成 taxId 在 NCBI 无对应记录），
不再"边搜边补 taxdump"，并且当命中表缺少 `staxids` 列时**直接报错**而不是静默丢命中。
无法定 taxId 的命中会写入 `unresolved_taxids.txt` + `unresolved_hits.tsv`。
另外，标准 BLAST `outfmt 6` 的 8 列表在 GTDB 模式下会被直接**拒收**：这种表没有
`staxids` 列，而它的 `qstart` 列会被误读成 taxId。

> **注意（上游既有行为，非 GTDB 分支引入）**：`search` 以**追加**方式写
> `<output>/<sample>.tsv`，因此对同一个 `-o` 目录重复运行会把每个蛋白写两遍，
> `analyze` 会把这些重复行都计入。重跑请换一个新的 `-o` 目录，或先把旧表改名保存。

## 分析

```sh
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank domain --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

- `--self-rank/--donor-rank/--tax-unirank` 接受 `domain`（自动映射到 taxdump 内实际存储的
  `superkingdom`）以及 `superkingdom` 别名。
- `--donor-label prefixed`（默认值即 `prefixed`，前缀指输出标签里 `g__` 这类官方前缀）输出官方 GTDB 标签（`g__Escherichia`）。
- `--donor-ncbi-map yes` 另存 `hgts/<sample>.donor_ncbi.tsv`，把 donor 映射回 NCBI taxId
  （来源：metadata 的 `ncbi_taxid` 列，r220 实测覆盖率 100%，完全离线）。
- `--bandwidth` / `predict.bandwidth` 可以是 `auto`、`grid`、`silverman`，也可以是 0.1–1.0 之间的数；
  写成 `bandwidth: 1`（YAML 整数）与 `--bandwidth 1` 现在得到同一个答案，
  越界值会被拒绝并写明"来源 + 可用区间"，而不是只甩一句 `Invalid bandwidth`。
- `--source` 的可用取值**按阶段不同**（由门禁逐对核对文档示例）：

  | 子命令 | 可用取值 | 默认 |
  | --- | --- | --- |
  | `database` | `ncbi\|gtdb\|hybrid` | `ncbi` |
  | `search` | `auto\|ncbi\|gtdb\|hybrid` | `auto` |
  | `analyze` | `auto\|ncbi\|gtdb\|hybrid` | `auto` |

  建库阶段没有 `auto`：那座库是用户**指定**要建哪一种，没有"从目录猜一个"的余地；
  而 search/analyze 面对的是既有库，才可以检测。写文档时不要把 `auto` 用到建库示例里。
- `--source` 在两条 analyze 路线上都按同一套判据生效：`ncbi` 不进入 GTDB 分支，
  `gtdb` 找不到 GTDB 标记就直接报错（点名 `gtdb/labels.tsv` 与出路），`hybrid` 找不到则打印
  WARNING 说明本次实际跑上游 NCBI 面，`auto` 才是"检测得到就用、检测不到就不用"。
- GTDB 分组的 `closesize` 默认上调为 200（古菌档 50），因为 GTDB 树远密于 NCBI；
  域由**每个蛋白最佳命中的域级祖先投票**自动判定（需严格多数，平局视为未判定），
  日志会打印 `GTDB profile: input domain = archaea -> using "gtdb_archaea"`。
  `--from-scores` 复用分数表时没有命中表，投票改用 `scores.tsv` 自己的 donor 列，
  因此两条路对同一个样本给出同一个档；同时会补一行 NOTE 说明分组阶段被跳过、
  这两个值本次不影响结果。检测不出域时日志写
  `input domain = not determined ... falling back to the default profile "gtdb"`，
  不会把"没检测出来"和"确实是细菌"混成同一句。显式给出的 `--close-size` 始终优先。
  可在 `config.yml` 的 `grouping.profiles` 调整。
- **参考库必须跨科以上**：实测用"肠杆菌科代表集"（793 基因组、单一科）配 E. coli 查询，
  4,311 蛋白在 `close_size=10/50/100/200` 下全部预测 0 条 HGT、无 distal 组——因为库里
  根本没有远缘类群。同理，查询基因组若已包含在库内，其最佳命中会被 self 组吸收。
  跨科/跨目抽样（如 `--rank-sample family`）才有统计意义。
- 实测（`db/r220_ar_family`，564 个古菌代表基因组，4,025 蛋白查询）：
  `close_size` 10/50/100 分别预测 371/371/132 条 HGT，说明该档位下古菌档 50 可用，
  而 200/400 需要接近全量代表集的库密度才有意义。

## GTDB-Tk 分类结果作为输入（P7）

GTDB-Tk **不是运行时依赖**（方案 N2）：本工具只读取它的结果表，把"用户基因组 → GTDB 谱系"
翻译成 GTDB 库的合成 taxId，供 `analyze --input-tax` 使用。

```sh
neoHGT gtdb gtdbtk --in BINs.classify.max_summary \
  --db "$NEOHGT_CACHE/db/r220_ar_family" \
  --out "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --report "$NEOHGT_CACHE/gtdbtk/gtdbtk_report.json" [--min-rank genus]

neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --input-tax "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --self-rank genus --donor-rank genus --donor-label prefixed
```

- 自动识别四种表形：v2 `classification`、v1 `classify_info`、de novo 的
  `gtdbtk.lineage` / `lineage`、以及"7 列等级"写法；v2 的 `#####` 横幅会被跳过，
  `#user_genome` 表头保留。
- **逐行回退**：v2 在有等级未指派时会把 `classification` 留空，程序自动改用
  `classify_info`，并把 NA 之后的段截断。
- taxId 取**最深的可匹配前缀**（例如只指派到属，就用属节点），报告 JSON 里有
  `matched_rank` 分布与每个未匹配基因组的**原因**。
- taxId 来源始终是目标库的 `gtdb/taxid.map.gz`；库里没有该谱系时判为 unmapped，
  不会凭空分配 ID。
- 实测（`db/r220_ar_family`，146,381 节点的 `taxid.map.gz`）：6/6 合成 gtdbtk 表行
  全部映射成功；查询基因组按属级指派 → `105000360`，随后 `analyze --input-tax`
  正常跑通端到端（`缓存测试/gtdbtk/`）。

## 谱系子集扩量到成员基因组（P7，可选）

GTDB 每个物种簇只发布一个代表基因组；本地图像只保存代表基因组的蛋白文件。
`--expand-members yes` 用 metadata 的 `gtdb_genome_representative` 列把"代表"展开成
"它代表的所有成员"（官方 `sp_clusters.tsv` 本地缺失，方案 §3.4 正是要求这样构造关系）：

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --expand-members yes --compile none \
  -o "$NEOHGT_CACHE/db/entero_members"
```

实测（r220）：793 个代表 → **94,393 个成员基因组**，全部写入
`gtdb/members.tsv.gz`（含 `representative` 与 `has_local_fasta` 两列）。
本副本**没有任何成员基因组的本地蛋白文件**（0/94,393），因此程序在抽取之前
就带着可操作提示中止（指出需要官方 `*_genomic_files_all` tar 或 `--seq-source
gtdb-tar|ena`），而不是悄悄产出一个更小的库。日志见 `logs/expand_members.log`。

## 系统发育距离辅助分组（P7）

本地发布包还带有**各 rank 的折叠树**（`tree_{bac120,ar53}/{x}_{rank}.tree`），
它们给出"到底有多少类群算近"的客观尺度，正是 `closesize` 想代理的东西：

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --top-k 20 --out dist.tsv --report dist.json
# 或直接指定树文件：--tree "$GTDB_DIR/tree_ar53/ar53_r220_genus.tree"
```

程序会用自研的轻量 Newick 解析器（无第三方依赖）算**patristic 距离**，
输出按距离排序的类群表与距离分布（`<= 0.01/0.05/0.1/0.2/0.3/0.5` 的类群数）。
`--gtdb-dir` 模式下会依次尝试该 rank 的两个域树，自动选中包含焦点类群的那一棵。

r220 实测（属级树）：

| 焦点 | 该树类群总数 | 中位距离 | ≤0.1 | ≤0.3 | ≤0.5 |
| --- | --- | --- | --- | --- | --- |
| `g__QDMS2`（古菌，本次查询基因组的属） | 1,846 | 3.00 | 0 | 1 | 27 |
| `g__Halorubrum`（古菌） | 1,846 | 3.03 | 1 | 10 | 37 |
| `g__Escherichia`（细菌） | 23,111 | 2.36 | 33 | 111 | 158 |

读法：古菌属级"真正近"的类群只有个位数到几十个（≤0.5 仅 27 个属），而细菌侧有 158 个属；
这与 §8 的 `close_size` 扫描完全一致 —— 古菌档必须比细菌/全量档小得多，
否则近邻组会吃掉全部最佳命中、distal 组为空。

**把树距离直接喂给 `analyze`**（phylogeny-assisted close group）：

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus --focus g__QDMS2 \
  --max-distance 0.5 --db "$DB" --out close.tsv --close-tax-out close.txt
neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$OUT" \
  --close-tax "$(cat close.txt)" --self-rank genus --donor-rank genus
```

`--db` 会用库的 `gtdb/labels.tsv` 把标签解析成合成 taxId，并只保留**本库基因组里
真实存在**的类群。另一层保护在 `analyze` 里：`--close-tax` 给出的 taxId 如果在命中表中
完全不可达（`analyze` 会先把分类树裁剪到命中涉及的范围），会被**丢弃并告警**而不是
让程序在读树时崩溃；若全部不可达则报错说明用法。实测（`db/r220_ar_family`）：
27 个 ≤0.5 的属中 21 个在该 564 基因归档无命中 → 被告警丢弃，close 组按剩余类群继续。
用解析出的 close 组在该档位重跑：预测从 371 条（9.2% 蛋白）提高到 849 条（21.1%），
而 53 个 ar53 保守标记基因仍全部得 0 分——即系统发育定义的 close 组没有把保守基因
误报为 HGT（同组数字见 README 的 GTDB 节）。

### 把供体清单映射回 NCBI taxId（`map-report`）

`analyze` 产出的 HGT 清单是逐行一个供体 taxId 或 GTDB 标签的文本；要离线地把它们
批量翻译回 NCBI taxId，用 `map-report`：它读库目录的 `gtdb/ncbi_map.tsv.gz`
（由 metadata 的 `ncbi_taxid` 列构建）与 `gtdb/labels.tsv`，逐行解析输入
（标签或合成 taxId 均可）：

```sh
neoHGT gtdb map-report --db <库目录> --in <HGT清单> -o <输出TSV>
```

`--to` 默认 `ncbi`（v1 只实现这一种目标体系）。输出 TSV 逐行给出映射结果，
并打印 rows / mapped / unmapped 三个计数；映射不到的行记为 unmapped，不算失败。

## 混合源 `--source hybrid`（实测）

`--source hybrid` 把 NCBI 分类树的顶部嫁接到 GTDB 树上：GTDB 的 `d__Bacteria` 与
`d__Archaea` 改挂到 NCBI 的 cellular root（taxId 131567），并按顶层 rank
（构造参数 `keep_ranks=('superkingdom',)`，别名表把当前 NCBI dump 写的 `domain`/
`kingdom` 一并认作顶层）保留若干真核与病毒顶层节点，其余 NCBI 节点不进树；
跨源共享蛋白的处置由 `--hybrid-lca-policy`（`root`/`domain`/`refuse`，方案 §5.7 第 3 项）
决定。NCBI 树的位置由 `--hybrid-ncbi-taxdump` 给出，缺失或没有 131567 都会明确报错
（提示改用 `--source gtdb`）。

> **这三档最初只是"被校验、被打印"，建库流程从未调用裁决函数（现已接入）**：抽取阶段遇到
> 重复蛋白 ID 时按"先遇到的 taxId 赢"写入 `taxon.map`，所以 `root` / `domain` / `refuse`
> 对真实构建没有任何区别，构建报告里那一行只是在复述用户输入（下表 IT-5 那条"三种 lca_policy
> 行为符合预期"验证的是合并后的树，不覆盖建库时的裁决，读它时请连同本段一起看）。现在的语义：
>
> - `root`：每个候选都保留，跨源共享蛋白的 LCA 落在 cellular root（RK7 记录的位移）；
> - `domain`：沿并树向上找到该候选自己的**域节点**（挂在 131567 下的 `d__Bacteria` /
>   `d__Archaea`，或 NCBI 侧被 graft 进来的顶层节点；新旧 dump 的 `domain`/`superkingdom`
>   拼写都认）并据此标注，多候选或跨域的蛋白记入 `x_source_proteins`，找不到域节点才落 root；
> - `refuse`：丢掉 taxId 不属于 GTDB 区间（`>= --gtdb-taxid-base`）的候选；若因此没有可留候选，
>   `taxon.map` 那一行写成 `na` 而不是猜一个。
>
> 裁决发生在抽取之后的一遍**流式替换**里（对 `taxon.map.gz` 的来源文件与 BLAST taxidmap），
> 内存只随"共享蛋白条数"增长。构建报告的 hybrid 一节现在会给真实条数：
> `proteins with more than one taxId candidate: N (reassigned M, left unassigned K)`——
> 纯 GTDB 序列的 hybrid 库里 N 通常是 0，这一行就老实写 0，也就是"本次构建没有可裁决的对象"。

用真实数据跑通的三条证据（GTDB 与 NCBI dump 均只读）：

| 门禁 | 做法 | 结果 |
| --- | --- | --- |
| IT-5 树合并 | `python tools/verify_hybrid.py --gtdb-dir $GTDB_DIR --ncbi-taxdump 缓存测试/input/ncbi_taxdump --out 缓存测试/reports/IT5_hybrid_r3.json` | PASS：纯 GTDB 树 146,380 节点，合并后 146,404（嫁接 25 个 NCBI 节点）；单根；两个域的父节点为 131567；细菌 vs 真核的 LCA 落在 cellular root（即 RK7 记录的语义位移）；三种 `lca_policy` 行为符合预期 |
| IT-5b 全链路 | 6 基因组 hybrid 库 → DIAMOND blastp → `neoHGT search` → `neoHGT analyze`（脚本 `临时文件/verify_hybrid6.sh`） | 与同 6 基因组的 `--source gtdb` 库**逐字节相同**：`hits.tsv` 相同、search 表数据行相同、`scores.tsv` 相同（42/42 命中带 staxid，HGT 判定同为 0）——往上嫁接不改变原核查询的结果 |
| IT-6e hybrid 断点续建 | `tools/verify_resume.py --fresh 临时文件/hybrid6_step1 --resumed 缓存测试/db/hybrid6`（第二次用全新空 `--tmpdir`） | PASS：7 个文件字节相同、6 个 gz 去注释后相同、5 个仅时间戳不同、0 不符、0 缺失 |

DIAMOND 的 rank 词表坑（实测发现）：`nodes.dmp` 的 rank 列会被 DIAMOND 校验，而
**当前** NCBI taxdump 用了新式 rank——`cellular root`(131567)、`acellular root`(10239
病毒)、`realm`、`clade`——`diamond makedb` 会直接 `Error: Invalid taxonomic rank:
cellular root` 退出。于是 `gtdb.write_taxdump()` 把不在 `gtdb.DIAMOND_RANKS` 词表里的
rank 落成 `no rank`（这正是 NCBI 在 2024 年改词表前对这些节点用的写法）。**只改文件**：
内存树、`gtdb/lineages.tsv`、`gtdb/labels.tsv` 保留真实 rank，`taxdump/nodes.dmp` 读出
来的树也仍按 taxId 与前缀标签工作。被压平的节点在 `meta.json` 的
`taxonomy.hybrid.ranks_flattened_for_diamond` 与 `gtdb/build_report.md` 的
"Hybrid graft" 小节逐个列出（真实数据上 `acellular root`×1、`cellular root`×1、
`clade`×1、`realm`×6）。GTDB 自己的 7 个 rank 经 `domain→superkingdom` 别名后全在词表内，
`gtdb.illegal_dmp_ranks()` 在纯 GTDB 树上恒为空，所以这条回退不会掩盖 GTDB 侧的问题
（`test_gtdb_ranks_are_never_flattened` 锁住；`test_modern_ncbi_ranks_survive_diamond`
锁住压平只发生在文件里）。

另一个实现约束：合并在**副本**上进行。早期实现直接改写被包装的 GTDB provider 缓存的树，
导致"纯 GTDB 树"里也混进了 NCBI 节点（谁先被问就决定它长什么样）；现在
`HybridProvider.load_taxonomy()` 先复制再嫁接，`test_merge_does_not_mutate_the_gtdb_tree`
锁住这一点。

仍未打通的部分：hybrid 库里**可搜的仍只有 GTDB 蛋白**（`--seq-source local`），因为
NCBI 蛋白序列的下载未获授权；嫁接进来的 NCBI 节点目前只服务树的语义。要让真核/病毒供体
真的能被命中（原论文 Asgard 场景所需），还得额外获取 NCBI 参考蛋白库——见 D12 §10.4
未决项。

## 资源预算（实测）

| 档位 | 基因组 | 蛋白 | db.faa | diamond db | 备注 |
| --- | --- | --- | --- | --- | --- |
| 6 基因组冒烟（IT-6 档，`fixtures/lists/e2e6.IDlist`） | 6 | 14,373 | 6.2 MB | 614 MB（base `1e8`） | 一次构建 81 s；续跑 99 s（跳过抽取） |
| 6 基因组冒烟（IT-2 档） | 6 | 27,667 | 12.3 MB | **620 MB**（base `1e8`）/ **5.6 GB**（base `1e9`） | makedb 峰值 RSS ≈ 3.8 GB；证明 `db.dmnd` 随 taxId **数值范围**膨胀 |
| `--rank-sample family --domain archaea` | **564** | **905,829** | **360 MB** | **985 MB** | 墙钟 **2 分 37 秒**，峰值 RSS 6.5 GB，可在 CI 运行 |
| `--gtdb-taxa f__Enterobacteriaceae --reps` | **793** | **3,247,172** | **1.41 GiB** | **2.07 GiB** | 墙钟 **4 分 50 秒**；`genomes.tsv` 793 行 = metadata 精确计数 |
| `--rank-sample genus`（细菌） | 23,112 | 6,000–8,000 万 | 28–37 GB | 42–110 GB | 按实测 466 B/蛋白外推 |
| 全量代表集 | 113,104 | 333,384,550 | **≈145 GiB**（外推） | ≈1.5–3× db.faa → 220–450 GB | 抽取按 793 档吞吐外推 ≈2–3 h + makedb 6–15 h；未执行 |

**体积基准**：实测 `db.faa` = **466 B/蛋白**（三档：455 / 417 / 466），容量预检
即按此值估算（峰值系数 3×）。方案给的"406 GB / 333M 蛋白"对应的是**全部 37.3 万个
每基因组文件**（含非代表集），因此直接用它外推会把代表集高估约 2.8 倍。

**跨卷读带宽实测**：单流 60 MB/s，8 并发 127 MB/s（NFS `/data2` → 本地 `/data`）。
低于方案设定的 200 MB/s 阈值 → **建议分档构建**（先 `archaea_family`，再 `genus`，最后
全量），并在构建后用 `--keep-faa compress` 释放空间。构建前会做容量预检
（空闲需 ≥ 3× 估算 `db.faa`），不足即中止，除非显式 `--force`。

**DIAMOND 的两个实测硬约束**（务必知悉）：

1. `db.dmnd` 体积随 taxId **数值范围**增长（而非节点数）：同一条 13 MB 库，
   `--gtdb-taxid-base 1e9` 得 5.7 GB（S1 spike 实测；IT-2 档用 12.3 MB 库复测得
   5.6 GB，见上表），`1e8` 得 620 MB。**因此本实现默认使用
   `100000000`**（仍比现网最大 NCBI taxId ≈3.4e6 高约 30 倍），并拒绝低于
   `10000000` 的取值。
2. `nodes.dmp` 的 rank 必须使用 NCBI 词表（GTDB 的 `domain` 会写成 `superkingdom`），
   且必须写满 13 列；否则 `diamond makedb` 直接失败。本实现已自动处理。

**BLAST 路线**：`makeblastdb` 对 local id 有 **50 字符上限**，而 GTDB 序列 ID 为
151–270 字符 → GTDB 模式下 `--compile blast|both` **直接报错退出**（方案 ADR-8），
只有显式加 `--force` 才会尝试执行。v1 只承诺 DIAMOND。

## 让 `--taxmap` 在完整规模下可用（P7，路线 B 复活）

方案原本**否决**路线 B（把 3.33 亿条"蛋白→taxId"读成 Python dict 需要数十 GB 内存）。
现在改为**磁盘锚点索引**：`gtdb.TaxidIndex` 只保留每 2048 行一个 (key, byte-offset) 锚点
（约为数据量的 0.2%），并把锚点表缓存在 `<map>.idx`；查询时二分定位 + 局部线性扫描，
**不再把映射载入内存**。

```sh
# 需要明文且按第一列字节序排序（即建库写出的顺序）
gunzip -k "$DB/taxon.map.gz"
neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  --taxmap "$DB/taxon.map" -o "$OUT" --tmpdir "$TMP"
```

- 命中表**没有** `staxids` 列、且序列 ID 也无法按路线 C 解析时，才由路线 B 兜底；
  解析顺序仍是 A → C → B → 报告 `unresolved_taxids.txt`。
- 传入 `.gz` 会被明确拒绝并提示先解压（gzip 无法按字节偏移随机访问）。
- 实测：10,000 行 / spacing 512 → 20 个锚点、索引文件 21 行；第二次打开直接复用缓存
  （<0.001 s）。测试：`TestTaxmapIndex` 5 项 + 工作流级
  `test_route_b_after_route_c`（无 `staxids`、路线 C 不认的 ID 由 `--taxmap` 补齐）。

## NCBI 路线的内存修复（P7，可选）

上游 `database.py` 的 `extract_genomes()` 会为**每一个唯一蛋白**在内存里保留一条记录，
并直接打印"可能会被 OS 杀掉"的警告；GTDB 分支从一开始就改成流式，本次把同等能力
回灌 NCBI 路线，作为**显式开关**（默认 `no`，默认路径的代码与行为完全不变）：

```sh
neoHGT database --source ncbi -c bacteria --stream-extraction yes -o "$DB" --tmpdir "$TMP"
```

- 去重改用 `--tmpdir` 内的 sqlite 索引，`genome.map.gz` / `taxon.map.gz` /
  BLAST `-taxid_map` / DIAMOND `prot.accession2taxid` 由**外部排序 + 流式分组**生成，
  峰值内存不再随蛋白总数线性增长；排序固定 `LC_ALL=C` 以保持与上游 `sorted()` 相同的字节序。
- 代价是 `--tmpdir` 需要额外空间（约每蛋白出现一次 30 字节）。
- 等价性由 `neoHGT/tests/test_database.py::TestStreamExtraction` 逐项断言：同一输入下
  `db.faa`、`genome.map.gz`、`taxon.map.gz`、每基因组蛋白数/残基数与上游实现**完全一致**；
  实测峰值 RSS 见 `benchmarks/stream_extraction_slope.txt`（两种规模 200k / 1M 蛋白）：

  | 实现 | 200k 蛋白 | 1M 蛋白 | 边际成本 | 外推到 333.4M 蛋白 | 抽取耗时(1M) | tmpdir(1M) |
  | --- | --- | --- | --- | --- | --- | --- |
  | 上游 | 365 MB | 1,071 MB | ≈882 MB/百万蛋白 | **≈287 GB** | 20.3 s | 0 |
  | 流式 | 193 MB | 304 MB | ≈138 MB/百万蛋白 | **≈45 GB** | 45.9 s | 104 MiB |

  即内存**降低约 6.4 倍**，代价是约 2.3 倍抽取时间与每百万蛋白约 100 MiB 临时盘。
- 需要 LCA 的共享蛋白在流式实现里**逐组**计算（组内只有几个 taxId），因此既省内存
  又与上游 `max(tids) / find_lca()` 的规则一致。

## KB 级测试夹具与离线端到端（§10.1 / IT-4）

`neoHGT/tests/data/GTDB/` 是一份**自成一个可搜索库目录**的 KB 级夹具（约 119 KiB），
由 `tools/make_fixtures.py` 从只读的本地 r220 数据生成：

```
GTDB/
├── protein/<6 个基因组的 .faa 摘录>     # 每个只留前 5 条蛋白（复制输入，不是软链）
├── metadata_bac120_sample.tsv / metadata_ar53_sample.tsv / taxonomy_sample.tsv
├── layout.txt / fixtures.json / meta.json
├── taxdump/{nodes.dmp,names.dmp}        # 转换器产物（32 个节点）
├── gtdb/{taxid.map.gz,labels.tsv,lineages.tsv,stem2taxid.tsv.gz,
│        genome2taxid.tsv.gz,ncbi_map.tsv.gz}
├── db.faa                                # 摘录合并（30 条蛋白）
└── search/sample.tsv                     # 7 列含 staxids 的命中样例（168 行）
```

> 这份清单不是"写过就算"：`tools/verify_fixture_provenance.py` 会重新打开只读来源，
> 逐条确认每条夹具序列确实存在于来源且残基一致、来源 sha256 与字节数没变、
> `protein/` 下没有未登记的孤儿文件、`db.faa` 的实际条数等于清单公布值、基因组数 ≤ 10。
> 报告见 `缓存测试/reports/R35_fixture_provenance.json`（实测：6 个来源全部读到、0 问题）。

重建（源数据只读；单次最多复制 10 个基因组，每个副本连 sha256 记入
`缓存测试/fixtures/manifest.tsv`；先 `--dry-run` 看待复制清单）：

```bash
. tools/env.sh
python tools/make_fixtures.py --dry-run
python tools/make_fixtures.py --limit 6
```

命中样例是**按谱系距离构造**的，不依赖本机 DIAMOND 版本，因此可跨环境逐字节复现：

- 多数蛋白以自身为最优命中，其后按"同科 → 同域 → 跨域"依次降分，用来演练 self/close 分组；
- 每个基因组末尾 2 条蛋白**没有 self 命中**，最优命中来自谱系距离最远的基因组，即人造 HGT。

`neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd` 用**真正的 CLI 解析器**
（`scripts/neoHGT` 的 `parse_args`）离线跑完整条链路，等价于：

```bash
neoHGT search  -i GTDB/protein/RS_GCF_000026225.1_…_g_Escherichia.faa \
               -o <out> -m precomp -s GTDB/search/sample.tsv \
               -t GTDB/taxdump --source gtdb --threads 1
neoHGT analyze -i <out> -o <out2> -t GTDB/taxdump --source gtdb \
               --input-tax <genome stem>:<species taxId> \
               --self-rank genus --donor-rank genus --donor-label prefixed \
               --close-size 1 --silhouette 0
```

断言的不是"跑通了"，而是**找回被埋的那 2 条**、donor 为 `g__Nanoarchaeum`、
`scores.tsv` 含 `distal` 列、`scatter.pdf` 与各 KDE/hist PDF 齐备、每条命中都有数字
taxId（DoD-3）、并且全程 `urlopen`/`run_command` 一次都没被调用（DoD-4）。
`--silhouette 0` 是因为 5 条蛋白的样本不足以支撑轮廓系数的判读，真实数据上不应关闭。
夹具在这两个测试里只是**只读来源**：测试先把它复制进临时目录再跑，末尾用整树 sha256
比对证明提交进仓库的夹具一个字节都没被写过。

## 零回归门禁（RT-1，贡献者必读）

GTDB 适配必须**不动 NCBI 路径的任何一个已有行为**，因此除了单测，还有一条机械可核对的门禁：

```sh
python tools/verify_regression.py \
    --baseline /path/to/pristine/neoHGT-main \
    --code     "$NEOHGT_CODE" \
    --out      "$NEOHGT_CACHE/regression/RT1_failure_set_diff.txt"
```

它做两件事：

1. **失败集合逐条比对**（不是比计数）：分别在两棵树跑 `neoHGT/tests`，按 test id
   求差；`newly_failing` 必须为空。上游本来就有 13 项失败（画图格式断言 + PATH 上
   没有外部二进制），该清单固化在 `ci/rt1_baseline_failures.txt`，CI 步骤直接比对它。
2. **删除行审计**：`neoHGT/{database,search,analyze,util}.py` 与 `config.yml` 相对
   基线必须 `deleted_lines = 0`（纯追加 + 方法重绑定）；唯一允许的改动是
   `scripts/neoHGT` 里 `module_names` 增加 `gtdb` 的那 1 行。

配套的 `tools/verify_ncbi_fixtures.py` 补上"夹具输出逐字节比对"那一半：在同一解释器、
同一随机种子下用两棵树各跑一次 NCBI 分析路径（`neoHGT/tests/data/Ecoli/search` →
silverman，再 `from_scores` → grid），逐文件 sha256 互比，并先各跑两遍证明树内可重复。
PDF 里 matplotlib 会写 `/CreationDate`、`/ModDate`，同一棵树两次运行的原始字节本来就不同，
因此 PDF 去掉这两个字段后比较，其余逐字节。注意：**仓库里已提交的
`Ecoli/analyze/scores.tsv` 在当前 scipy/scikit-learn 下连原始上游代码都复现不出来**
（两棵树产出相同、但与夹具不同，上游 `test___call__` 也因此在基线里同样失败），
所以这条门比的是"两棵树是否一致"，不是"是否等于那份旧夹具"。

还有一条"名实相符"的门禁：`tools/verify_plan_traceability.py` 从方案文档里抽出所有
`UT-*/IT-*/AT-*/RT-*` 验证条目编号，要求每个编号要么被某个测试文件署名引用，
要么映射到一个**确实存在**的证据产物；否则退出码 1。产物是
`缓存测试/reports/plan_traceability.md` 这张对照表。这条门是有效的：UT-0.3、
UT-2.1/UT-2.x、UT-3.3、UT-3.4、UT-4.x、UT-5.x 都曾是"方案里写着、实际没有对应测试"，
由它扫出并补齐（AT-4 更早由人工翻方案发现，纳入同一张表核对）；把方案里任一编号改成
不存在的值即可看到它报 FAIL。

如果 pytest 没有输出汇总行（解释器缺失、收集错误等），工具**拒绝比较**而不是报告通过——
空集合互相"相等"正是最容易伪装的假绿。CI 里同理：断言基线文件非空、断言汇总行存在、
断言解析出的 `FAILED` 行非空，任何一条不满足都判失败。

## 已知限制

### 分组参数的实测注意事项（三个样本）

`缓存测试/reports/S10_4_3_sensitivity.md` 是在 halo（古菌）、bathy（不同门的古菌）、
entero（细菌，单科库）三个样本上重扫 `close_size` / `input_cov` / `bandwidth` / `noise` 的结果。
三条会改变结论的注意事项：

1. **`--close-size` 是下限，不是配额。** 默认按域取档（细菌 `grouping.profiles.gtdb` = 200、古菌 `gtdb_archaea` = 50），在树上取不出那么多类群时会被截断（bathy 在它的默认档只有 20 个 close 类群），所以"用了默认值"不等于"close 组有 200/50 个"。
   日志里 `Close group has N taxa.` 那一行才是实际生效值，跨样本比较必须看它。
2. **`--noise` 的方向不保证单调。** 配置注释写着"越大越严格"，但实测两个样本方向相反
   （halo 25→75：374→325 变严；bathy 25→75：285→421 变松）。它改变的是谷-峰距离里被排除的
   区间端点，在双峰或长尾分布上方向会翻。调这个参数必须看自己样本的曲线，不能引用注释。
3. **`--input-cov` 的作用是输入分类，不是预测条数。** 命中分类干净时（最优命中 100% 落在同一
   类群）50–95 全平坦；只有覆盖率处在阈值附近（细菌样本实测 89.47%）才在 cov=95 时把输入分类
   从种退到属。因此它属于"样本本身分类是否可信"的开关，`--self-rank` 决定这个退化是否影响结果。
4. **`donor=0` 不是"没有供体"。** 它表示远缘命中的 LCA 粗于 `--donor-rank` 所请求的秩：
   E. coli 对 4,896 基因组跨科库的实测是，请求属级时 10.4% 可标注、科级 13.1%、域级 100%，而
   那些 `0` 的 LCA 节点本身有 591 个落在 superkingdom、65 个 class、21 个 phylum。`--distal-top`
   在这个数据集上 10–100 全平坦，不是可用的调节旋钮；要更高解析度就放宽 `--donor-rank`。
5. **`--from-scores` 与整跑必须同结果，但曾经不同。** GTDB 的 profile 默认值与
   `domain → superkingdom` 秩翻译原先只做在 `read_input()` 里，而 `--from-scores` 跳过该步骤，
   于是 `--donor-rank domain` 在复用路径上一条都标不出来。现已抽成 `apply_gtdb_grouping()` 由两条
   路线共用（实测修好后两路完全一致：1235 条预测、100% 标注），并由
   `neoHGT/tests/test_gtdb_pipeline.py::TestFixtureEndToEnd` 钉住。
6. **同一座库、同一个参数，四个细菌样本的 `close_size` 曲线并不一致（实测）。**
   对 4,896 基因组"每科一个代表"的库：*Escherichia* 10/50/100/200 恒 771、Alicyclobacillaceae
   恒 753、f__ATN3 恒 714，而 **Arcobacteraceae 在 10/50 给 436 条、到 100/200 直接归零**——
   它的 self 组只有 2 个类群（该属在这座库里只有一个别的代表），`close_size≥100` 时推断出的近邻组
   一口气涨到整个门（13,208 个类群），把全部命中吞掉，distal 组为空，
   日志给出 `WARNING: No hit is assigned to distal group. Cannot predict HGTs.`。
   所以"细菌档 200 可用"是**逐样本**的性质，不是库的性质；数据取自
   `缓存测试/reports/R26_bac_close_size_curves.json`（纯预测档，未重跑 DIAMOND）。

顺带确认：`--from-scores` 复用已落盘的 `scores.tsv` 与整跑**在预测集合上同值**（bathy 两档均为
724→367），所以做阈值类参数（`bandwidth`/`noise`/`fixed`）的扫描可以放心走复用路径，每次省掉整轮打分。
但"同值"只对被调用的蛋白集合与排名成立，浮点列会因 `scores.tsv` 的 6 位有效数字落盘差到最后一位
（见本小节第 5 条），比较时不要断言逐字节相等。

- 无真核/病毒：GTDB 仅 `d__Bacteria` 与 `d__Archaea`（实测）。
- hybrid 模式需要 NCBI taxdump：本地 GTDB 副本不含该文件，本次经使用方单独授权后才下载
  `input/taxdump.tar.gz`。**混合树已用真实 NCBI taxdump + 真实 GTDB r220 联调通过**：
  IT-5（树合并与 LCA 语义）、IT-5b（6 基因组建库 → DIAMOND 搜索 → analyze 全链路）、
  IT-6e（hybrid 断点续建）。仍未获授权的是 hybrid 的**NCBI 蛋白序列**通道
  （`--seq-source ncbi`），因此 hybrid 库里可搜的仍只有 GTDB 蛋白——嫁接进来的 NCBI
  分类节点目前只服务树的语义，不贡献命中；要真做出原论文 Asgard 场景所需的"真核/病毒
  供体可被命中"，需要额外获取 NCBI 参考蛋白库（见 D12 未决项）。
- 远端序列通道（`--seq-source ncbi|gtdb-tar|ena`）只导出抓取计划
  `gtdb/fetch_plan.tsv` 供人工执行——本版本没有下载执行器，不会发起任何联网下载；
  `--allow-network` 默认按 no 处理（不是联网开关），`--allow-network yes` 只是使用这些通道的前置确认——yes 不等于会联网，下载执行器不存在。
- 文件名内嵌谱系止于属级，种级一律以 metadata/taxonomy 为准。
- 本地副本经第三方脚本改名、重写序列 ID，并存在与 r207 合并的痕迹；无官方 MD5SUM，
  故本工具自建输入 sha256 指纹并强制自洽核对。
- `--keep-faa drop` 当前**只告警不删除**（删除需先备份，符合项目安全规则）。
- GTDB taxId 是合成 ID，不可与 NCBI taxId 混用于同一棵历史结果树；跨 release 依赖
  `gtdb/taxid.map.gz` 保持 ID 稳定——这一点已在本地 r214/r220 之间双向实测（0 改动 / 0 冲突 /
  0 复用，见"单独构建 taxdump"一节末与 `缓存测试/reports/R24_lineage_only_and_release_upgrade_SOP.json`）。
- 官方 release 的**扁平布局**（谱系表直接放在数据根目录）现在会被识别，`gtdb taxdump` 也能只凭
  谱系表工作；但完整建库仍需 metadata 表，而本地另外两份副本（r214、r226）都不带 metadata，
  因此"从官方数据集一次跑到底"仍未实测过。列映射本身是照官方 schema 校准的：`REQUIRED_COLS`
  的 11 列在真实 `bac120_metadata_r220.tsv` 表头里 11/11 命中。
- r226 那份本地副本的谱系串被第三方用 `_` 粘连并截断到属级，属于**不合格输入**：深度核对会在
  任何写入之前拒绝它（退出码 1），需要改用官方原始的 `*_taxonomy_r226.tsv`。
