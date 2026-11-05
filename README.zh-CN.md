[English](README.md)(主)| **简体中文**(与英文版内容对照)


# neoHGT2GTDB: 🧬 基于 GTDB 参考框架的水平基因转移

> **这不是官方的 neoHGT / neoHGT2。** neoHGT2GTDB 是 [neoHGT2](https://github.com/cgneo/neoHGT) 的一个独立、非官方的适配版本，在本仓库中维护。它未经上游作者或 [Qiyun lab](https://github.com/qiyunlab) 的认可，与其没有任何隶属关系，也不受其支持。官方项目仍位于 [cgneo/neoHGT](https://github.com/cgneo/neoHGT)；此处不使用任何官方 logo 或品牌标识。

**neoHGT2GTDB** 是对 [neoHGT2](https://github.com/cgneo/neoHGT) 的一次深入、以 GTDB 为核心的适配。它原样保留了上游项目——NCBI 路线、多进程数据库下载器、搜索与分析流水线——并在此基础上添加了一条完整的 GTDB 路线：由本地只读的 [GTDB](https://gtdb.ecogenomic.org/) 发行版（已在 R09-RS220 上验证）构建的参考数据库、记录在只追加清单中的合成 taxId、GTDB 风格的阶元与供体标签、可选的 NCBI+GTDB 混合树，以及为上述一切而设的机器校验文档/测试门禁。

命令行入口点仍是 `neoHGT`，继承自上游，因此既有的 neoHGT/neoHGT2 脚本可以继续工作；GTDB 方面的新增内容位于 `gtdb` 子命令实用工具以及新的 `--source gtdb|--hybrid` 数据库选项中（参见 [GTDB 数据库（实验性）](#gtdb-数据库实验性)）。

**命名，在此一次说清：** *仓库* 和 *pip 发行包* 的名字是 **neoHGT2GTDB**；*导入包* 和 *命令* 仍为 `neoHGT`。这是一个深思熟虑的兼容性选择：每一个文档示例、既有的 neoHGT/neoHGT2 脚本以及机器校验的文档门禁都固定使用上游的拼写，重命名可执行文件只会破坏这一切而不会带来任何好处。因此 `pip install .` 会安装一个名为 `neoHGT2GTDB` 的发行包，它提供 `neoHGT` 命令和 `neoHGT` 包。

项目沿革：[HGTector](https://github.com/qiyunlab/HGTector) → neoHGT (v1) → [neoHGT2](https://github.com/cgneo/neoHGT) → **neoHGT2GTDB**（本仓库）。如果你有兴趣复现[论文](https://www.nature.com/articles/s41467-023-36487-z)中的结果，前一版本 neoHGT1 仍可在[这里](https://github.com/cgneo/neoHGT-v1)找到。继承自 neoHGT2 的 NCBI 路线中的 bug 最好报告到[上游的 issue 跟踪器](https://github.com/cgneo/neoHGT/issues)；GTDB 适配中的 bug 则应报告到本仓库的 issue 跟踪器，欢迎合作者加入。

## 目录
1. [关于](#关于)
2. [安装](#安装)
3. [数据库](#数据库)
4. [编译](#编译)
5. [GTDB 数据库（实验性）](#gtdb-数据库实验性)
6. [搜索](#搜索)
7. [分析](#分析)
8. [卸载](#卸载)
9. [许可证](#许可证)

## 关于
**neoHGT2GTDB** 为感兴趣的物种预测水平/侧向基因转移（HGT）。与它所源自的上游 [neoHGT](https://github.com/cgneo/neoHGT) 项目一样，它基于 Qiyun lab 开发的 [HGTector](https://github.com/qiyunlab/HGTector)——上游修复了 HGTector 的已知 bug，用多进程并行化了下载过程，并提升了整体稳定性。neoHGT2GTDB 自身添加的是 GTDB 路线：它在 GTDB 标准化、阶元完备的原核生物分类体系上重建参考框架（并带有可选的混合模式，把 NCBI 树的顶端嫁接到 GTDB 树上），这可度量地改善了 self/close/distal 分组的统计行为，并对该行为进行端到端的文档化和机器校验。

  

## 安装
本项目完全以 Python3 编写。<br>
推荐安装方法：

第 1/3 步：

```
conda create -n neoHGT -c conda-forge python=3 pyyaml pandas matplotlib scikit-learn
```

第 2/3 步：

接下来，你需要手动安装 diamond，方法是访问 [conda 官方网页](https://anaconda.org/bioconda/diamond/files?version=0.9.26)。

你会看到一个文件列表，请查看文件名的开头，它描述了操作系统。下载与你的系统相对应的那个文件后，你可以运行

```
conda activate neoHGT
```
```
conda install <diamond___.tar.bz2>
```

然后接上你刚下载的文件（把文件拖到终端里），即可安装 diamond。

第 3/3 步：

运行以下命令：

```
conda activate neoHGT
```
```
pip install .
```

请在本（neoHGT2GTDB）仓库的克隆中运行此命令。若要改为安装上游项目，请使用 `pip install git+https://github.com/cgneo/neoHGT.git`。


## 数据库

与 HGTector 相比，本项目继承的 neoHGT2 下载器不仅修复了 HGTector 的一些已知 bug，还使用多进程实现并行下载，为你节省数据库下载过程的时间。

neoHGT2GTDB 保留了上游的 ```server.py```，它无需图形用户界面即可下载并安装数据库。

```server.py``` 是**仅支持 NCBI** 的：它获取预构建的 `hgtdb_20230102`
数据库，对 GTDB 一无所知。GTDB（或 `--source hybrid`）数据库要
通过命令行构建——参见下文的
[GTDB 数据库（实验性）](#gtdb-数据库实验性)。

但推荐的方式是手动下载：

## Option1（手动）
### （最全面、最新但耗时）：

首先，我们需要 ```taxonmap```，获取该文件有两种选择：

下面的链接会得到一个更小、更紧凑的 ```taxdump```（解压后约 0.4GB）。你可以把链接粘贴到浏览器里开始下载：
```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
```

下面的链接会得到最新、更全面的 ```taxdump```（解压后约 1.8GB）：

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z
```

___

接下来，我们需要 ```prot.accession2taxid```，它有两种选择：

下面的链接会得到一个更小、更紧凑的 ```prot.accession2taxid```（约 7GB）。你可以把链接粘贴到浏览器里开始下载：

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

下面的链接会得到最新、更全面的 ```prot.accession2taxid```（约 12GB）：

```
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
```

## Option2（手动）：
### （省时但不够全面、不够新）

在[这个 OneDrive 链接](https://arizonastateu-my.sharepoint.com/personal/qzhu44_asurite_asu_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fqzhu44%5Fasurite%5Fasu%5Fedu%2FDocuments%2FPublic%2FHGTector)上有一个预构建好的数据库，你只需下载最新的 ```hgtdb___.tar.xz``` 文件，即可跳过数据库编译步骤，直接跳到 [编译](#编译) 步骤。

___

## 数据库下载

在获得 ```taxonmap``` 和 ```prot.accession2taxid``` 之后，我们需要运行 ```neoHGT database``` 命令，下载约 291,000 个基因组文件来构建 diamond 数据库（~90GB），具体取决于你选择了多少个类别。

构建数据库的命令是：

```
neoHGT database -c bacteria,archaea,viral -o <FOLDER_NAME>
```

其中 ```-c``` 或等价的 ```-cats``` 表示类别。

下表解释了可供你构建自定义数据库的命令和选项。

| 命令 | 选项 |
|---------|---------|
| -c (-cats)      | archaea, bacteria, fungi, invertebrate, metagenomes, mitochondrion, plant, plasmid, plastid, protozoa, unknown, vertebrate_mammalian, vertebrate_other，以及 viral（源自 [NCBI RefSeq 数据库的文件夹名](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)） |
| -o      | 输出数据库目录 |

___

已知问题：

如果看到这条错误消息，不要慌张：
```
File "/Users/___/miniconda3/envs/neoHGT/lib/python3.11/site-packages/neoHGT/database.py", line 351, in get_categories
    raise ValueError(
ValueError: "___" is not a valid RefSeq genome category
```

出现这个错误，要么是你输入的类别名称有误（例如 NCBI 服务器上并不存在 ```virus``` 文件夹，你应该输入 ```viral```），要么是 NCBI 服务器暂时切断了与你的连接（拒绝了你的连接）。程序本身没有任何问题。

遇到这个错误时，你能做的最好的事就是在终端输入以下命令来测试连接：

```
rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```

如果你看到类似 ```rsync: failed to connect to ftp.ncbi.nlm.nih.gov: Connection refused (61)
rsync error: error in socket IO (code 10)``` 的输出，这就证实确实存在网络问题（很可能出在 NCBI 服务器一侧）。

否则，如果前一条命令返回的列表以 ```drwxr-sr-x``` 开头，后面跟着类似 ```lrwxrwxrwx``` 的条目，就说明连接已恢复，你可以放心地再次执行程序。

关于更高级的选项（例如排除特定 taxid），请参考 HGTector 的[数据库](https://github.com/qiyunlab/HGTector/blob/master/doc/database.md)文档。

看到程序以下面的输出结束时：
```
✅ Done
Combined protein sequences written to db.faa.
```
即表示你可以继续进行手动编译。

___

## 编译

编译只使用 DIAMOND，因为它经过优化、速度更快。

### 如果你按照 <b>Option1</b> 的说明操作：

```
cd <DATABASE_DIRECTORY>
```

```
mkdir diamond
```

然后，仔细地用你下载并解压好的文件替换 ```<‼️>``` 处的条目。确保这些文件都位于你的 ```<DATABASE_DIRECTORY>``` 的根文件夹下。下面为了清晰起见，命令被拆分成多行，但在实际操作中你应该保持为一行，不要按 ```return``` 键。

```
diamond makedb --threads 64 --in db.faa
--taxonmap <‼️prot.accession2taxid.FULL>
--taxonnodes <‼️taxdump>/nodes.dmp
--taxonnames <‼️taxdump>/names.dmp
--db diamond/db
```

### 如果你按照 <b>Option2</b> 的说明操作：

第 1/2 步：

```
cd <DATABASE_DIRECTORY>
```

```
echo $'accession.version\ttaxid' | cat - <(zcat taxon.map.gz) > prot.accession2taxid.FULL
```
在这一步，你可能会遇到一个奇怪的错误，例如 ```zcat: can't stat: taxon.map.gz (taxon.map.gz.Z): No such file or directory```。如果看到这个错误，只需把 ```taxon.map.gz``` 文件重命名为 ```taxon.map.gz.Z``` 即可。

第 2/2 步：

确保你仍在 ```<DATABASE_DIRECTORY>``` 中，然后输入以下命令：

```
mkdir diamond
```
```
diamond makedb --threads 64 --in db.faa --taxonmap prot.accession2taxid.FULL --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp --db diamond/db
```

## GTDB 数据库（实验性）

这是 neoHGT2GTDB 的核心，也是相对于 neoHGT2 的深层改进：
除了继承自上游的 NCBI 路线之外，它还可以直接从**本地的只读 GTDB 发行版**（已在
R09-RS220 上验证）构建参考数据库。什么都不会被下载：分类体系由 GTDB 自带的
`*_taxonomy_r220.tsv` 加上
`*_metadata_r220.tsv` 转换而来，蛋白序列从每个基因组的
`.faa` 文件流式读取，并且每个合成 taxId 都记录在一份只追加清单中，
使 id 在多次重建之间保持稳定（`database` 和
`neoHGT gtdb taxdump` 在清单存在时会复用输出目录的清单，或
`--gtdb-id-manifest` 指定的清单）。

无需任何工作区设置：输出可以放在任何地方，除非
设置了 `$NEOHGT_ROOT`，那样写入就必须留在其中。每种配置中唯一一条硬性
规则是：只读的 GTDB 输入数据集绝不能被写入。远程序列通道
（`--seq-source ncbi|gtdb-tar|ena`）
在本版本中**没有下载执行器**：它们会写出
`gtdb/fetch_plan.tsv` 供手动执行，然后拒绝继续，而不是
悄悄构建一个空数据库。

这些实用工具位于 `neoHGT gtdb <subcommand>` 之下
（`probe, taxdump, map-report, gtdbtk, reconcile, tree-dist`），它与
`python -m neoHGT.gtdb <subcommand>` 是同一接口；较旧的扁平形式
`neoHGT gtdb --action probe|taxdump|map-report` 仍可解析。

```sh
# 0) read-only pre-flight: layout detection + 19 consistency checks + sha256 fingerprints
neoHGT gtdb probe --dir /path/to/GTDB_r220 --out "$NEOHGT_CACHE/logs/layout_probe.json"
# 0b) same checks, but write the accession-level difference lists to disk
neoHGT gtdb reconcile --dir /path/to/GTDB_r220 --out "$NEOHGT_CACHE/reconcile" \
  --report "$NEOHGT_CACHE/reports/reconcile.json"

# 1) development tier: one representative genome per family (564 archaeal genomes, ~2.5 min)
neoHGT database --source gtdb --gtdb-dir /path/to/GTDB_r220 \
  --rank-sample family --gtdb-taxa d__Archaea --seq-source local \
  --compile diamond --threads 16 -o "$NEOHGT_CACHE/db/ar_family"

# 2) search (the DIAMOND db carries the GTDB taxIds, so staxids come out directly)
diamond blastp --query sample.faa --db "$NEOHGT_CACHE/db/ar_family/diamond/db" \
  --threads 16 --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive -o hits.tsv
neoHGT search -i sample.faa -m precomp -s hits.tsv \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" -o "$NEOHGT_CACHE/search/sample"

# 3) analyse with GTDB-flavoured ranks and donors
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank genus --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

不兼容的 CLI 缩写（argparse 前缀匹配）：新的 GTDB 选项
与上游长选项共享前缀，因此 `database` 上的缩写
`--a --g --ge --gen --l --rep --ty --typ --type`、`search` 上的
`--r --re --s` 和 `analyze` 上的 `--donor-n` 变得有歧义，
现在会报 “ambiguous option” 错误而失败。请把长选项完整写出；
这是对既有 NCBI 脚本唯一的一处破坏性变更。

关于与 NCBI 路线不同行为的说明：

- 远程分类信息获取被**强制关闭**（`--fetch-enable no`），因为
  合成 GTDB taxId 在 NCBI 中没有记录。在 GTDB 路线上，缺少
  `staxids` 列的命中表会被大声拒绝，而不是悄悄丢弃
  命中（`--require-staxids` 在这里默认为 yes）；在 NCBI 路线上
  则保留上游行为。8 列的表（标准 BLAST
  `outfmt 6`）在 GTDB 模式下会被拒绝：它的 `qstart` 列会被读作
  一个 taxId。
- `nodes.dmp` 以 NCBI 阶元词汇写出（`domain` 被存为
  `superkingdom`），共 13 列，这是 `diamond makedb` 所要求的；在 CLI 上
  `domain` 和 `superkingdom` 均可接受并自动映射。
- 供体以官方 GTDB 标签（`g__Halorubrum`）报告，
  `--donor-ncbi-map yes` 还会额外写出 `hgts/<sample>.donor_ncbi.tsv`，
  离线地依据 GTDB 元数据的 `ncbi_taxid`
  列把每个供体映射回它的 NCBI taxId。
- `neoHGT gtdb reconcile --dir <release> --out <dir>` 对元数据、GTDB
  分类和本地 `.faa` 索引做差分，并写出五份 accession 级差异列表
  （`missing_faa`、`extra_faa`、`nonrep_faa`、
  `metadata_not_in_taxonomy`、`taxonomy_not_in_metadata`），即使它们为空也照写，
  从而留下审计痕迹，而不是只有一个通过/失败计数；任何差异都会以非零状态退出。
  在这份 r220 副本上，全部五个列表均为空
  （596,859 个基因组 / 113,104 个代表基因组 / 113,104 个 FASTA，24 秒，1.7 GB）。
- `gtdb/` 下的每个派生表（以及 `lineages.txt`、reconcile 列表、
  fetch 计划、成员扩展）都以一行溯源信息开头：
  `# neoHGT version=…; gtdb_release=…; layout=…; generated=…`，因此任何产物
  都能追溯到产生它的代码、发行版和时刻。
  `taxon.map.gz` 和 `genome.map.gz` 则有意保持无注释：DIAMOND、
  BLAST 和 `read_prot2taxid` 要原样使用它们。
- GTDB-Tk 的输出可以用作输入分类，
  而不必把 gtdbtk 变成依赖：`neoHGT gtdb gtdbtk --in BINs.classify.max_summary --db <db>
  --out input_tax.tsv` 会把 GTDB-Tk 谱系（v2 的 `classification`、v1 的
  `classify_info`、从头计算的 `gtdbtk.lineage`，或七个阶元列）
  转换为该数据库的合成 taxId，匹配最深的已分配阶元，
  其结果可直接送入 `neoHGT analyze --input-tax`。
- GTDB 会在不同发行版之间更改高阶分类群的名称（`p__Proteobacteria` →
  `p__Pseudomonadota`），而官方的 `synonyms.*.tsv` 并不包含在发行版
  副本中，因此 `neoHGT analyze --taxon-synonyms aliases.tsv`（`alias<TAB>current`）
  可以把你一直沿用的旧名映射到当前名称，作用于
  `--self-tax/--close-tax/--tax-include/--tax-exclude/--input-tax`，并打印
  哪个已弃用的名称被改写了（计划风险 RK11）。
- `--keep-placeholder-names`（默认 `yes`，来自
  `taxonomy.profiles.gtdb.keep_placeholder`）之所以存在，是因为有 1,658 个
  r220 分类群使用源自 accession 的代码作为官方标签（`2-12-FULL-45-22`、`g__DUMJ01`、
  `s__MK2026 sp010841923`）。正是出于这个原因，GTDB 配置默认关闭
  `--tax-capital`/`--tax-latin`；如果你开启其中之一，这些代码仍会被保留，
  除非你传入 `--keep-placeholder-names no`。
- 在 GTDB 下，`-m precomp` 指向一个包含逐样本表的**目录**会被拒绝：
  上游的循环给每一批都喂第一张表（其自身代码中自带的
  `# todo: general fix`），因此后面的每个样本都会悄悄返回零命中——
  在使用合成 GTDB taxId 时，这看起来像“没有同源序列”而不是“输入
  文件有误”（计划附录 F6，已登记于 `docs/既有缺陷登记.md`）。
- `--expand-members yes` 会把一份代表基因组的选择扩展为
  它所代表的物种聚类的成员基因组（由于本地缺少 `sp_clusters.tsv`，
  关系取自元数据列 `gtdb_genome_representative`），
  并把完整列表写入 `gtdb/members.tsv.gz`。在这份 r220 副本上，793 个
  Enterobacteriaceae 代表基因组会扩展出 94,393 个成员，而它们都不附带
  蛋白文件，因此构建会带着一条可操作的提示停下来，而不是
  悄悄交付一个更小的数据库。
- `neoHGT gtdb tree-dist --gtdb-dir <release> --rank genus --focus g__X`
  会依据该发行版的逐阶元树报告 patristic 距离，客观地度量有多少分类群是
  “接近”的。古菌的属级树中，只有 1,846 个属中的 27 个落在查询属
  的 0.5 次替换以内（细菌方面：Escherichia 为 23,111 个属中的
  158 个）——这正是 `gtdb_archaea` 层级所需的 `closesize` 要比细菌或全量层级小得多的原因。加上 `--db <db>
  --close-tax-out close.txt` 可以把这些标签解析为该数据库的合成
  taxId，从而可直接传给 `neoHGT analyze --close-tax`，即一个
  由系统发育定义的 close 组；在古菌层级上，这把预测
  从 371（占蛋白的 9.2%）提高到 849（21.1%），
  同时 53 个保守的 ar53 标记基因中
  仍有 0 个被判为 HGT。
- `--stream-extraction yes`（仅限 NCBI 路线；默认关闭，行为不变）
  用一种基于磁盘的等价实现替换上游的蛋白提取：
  实测每百万蛋白 882 -> 138 MB，即对 333M 个蛋白预计
  ~287 GB -> ~45 GB，提取时间为 2.3x，每百万蛋白约 ~100 MiB
  临时磁盘，且 `db.faa` / `genome.map.gz` / `taxon.map.gz` 逐字节相同。
- 重新运行完全相同的命令会复用 `db.faa` 和已记录的
  `gtdb/extraction_state.json`，因此一次中断的数小时提取
  不会重复进行；重新生成的数据库除时间戳外与一次成功的构建
  逐字节相同。
- `neoHGT search` 会**追加**写入 `<output>/<sample>.tsv`（上游行为，
  此处未改动）：向同一个 `-o` 目录运行两次会把每个蛋白写两遍，
  而 `analyze` 会对两份副本都打分。每次运行请使用全新的 `-o`
  目录，或者先把上一次的表移开。
- `db.dmnd` 随 taxId 的数值*范围*增长，因此 id 基数默认为
  `100000000`（`--gtdb-taxid-base`，最小值 `10000000`）。
- `--source hybrid` 把 NCBI 树的顶端嫁接到 GTDB 树上，使真核和病毒分类群
  能与 GTDB 原核生物共用一棵树：两个 GTDB 域被重新挂到 cellular root
  （taxId 131567）之下，
  `--hybrid-ncbi-taxdump` 指向 NCBI 的 `nodes.dmp`/`names.dmp`（必需），
  而 `--hybrid-lca-policy root|domain|refuse` 决定两个来源都认领的蛋白
  会如何处置。当前的 NCBI dump 把最高阶元写作
  `cellular root`、`acellular root`、`realm` 和 `clade`，而 DIAMOND 在 `makedb` 处
  会拒绝其中每一个，因此 `nodes.dmp` 会把词汇表之外的任何阶元折叠为 `no rank`，
  而合并后的树、`gtdb/lineages.tsv` 和
  `gtdb/labels.tsv` 保留真实阶元；被折叠的内容记录在
  `meta.json` -> `taxonomy.hybrid` 以及 `gtdb/build_report.md` 的 “Hybrid graft”
  一节中。已在真实数据上验证：合并本身
  （`tools/verify_hybrid.py`：25 个 NCBI 节点嫁接到 146,380 个 GTDB 节点上，
  单一根，跨来源 LCA 位于 cellular root），以及一个 6 基因组混合数据库——
  其 DIAMOND 命中、搜索表和 `scores.tsv` 与用 `--source gtdb` 构建的相同 6 个基因组
  逐字节一致。可搜索的序列仍然只有 GTDB 蛋白，
  因为在本环境中下载 NCBI *蛋白*序列
  未获授权。
- BLAST 路线**无法**用于 GTDB id：`makeblastdb` 把本地 id 限制在
  50 个字符以内，而 GTDB 序列 id 长达 151-270 个字符。
  仅支持 DIAMOND。

实测开销（R09-RS220，全集含 113,104 个代表基因组 / 333,384,550
个蛋白）：

| 层级 | 基因组数 | 蛋白数 | db.faa | db.dmnd | 墙钟时间 |
| --- | --- | --- | --- | --- | --- |
| 6 个基因组（冒烟测试；`缓存测试/db/it6` 一类的构建） | 6 | 14,373 | 6.2 MB | 614 MB | ≈1 分钟（81 秒；续跑 99 秒） |
| 科级抽样，古菌 | 564 | 905,829 | 360 MB | 985 MB | 2 分 37 秒 |
| 一个科，`--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 分 50 秒 |
| 属级抽样，细菌 | 23,112 | ~60-80 M | 28-37 GB | 42-110 GB | 外推 |
| 全部代表基因组集合 | 113,104 | 333,384,550 | ~145 GiB | ~1.5-3x db.faa | ~2-3 小时提取 + 6-15 小时 makedb（外推） |

实测：`db.faa` 每个蛋白占用 **466 字节**（三个层级：455 / 417 / 466），
磁盘预检正是按此假设，并为 `makedb` 留出 3x 的峰值系数。

中文文档（包括工作区约定和完整
参数列表）位于 [`docs/GTDB-数据库.md`](docs/GTDB-数据库.md)。

## 搜索

下面这种方法被认为是一种更“安全”的方式，因为你将使用 ```diamond``` 或 ```blastp``` 手动执行搜索，这可以降低因第三方库支持（例如 diamond 方面）而出现 bug 的概率。

第 1/3 步：

```
conda activate neoHGT
```

第 2/3 步：

（下面为了清晰起见，命令被拆分成多行，但在实际操作中你应该保持为一行，不要按 ```return``` 键）

```
diamond blastp
--query <input_file.faa>
--db <diamond_db>
--outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids
> <output_directory>/result.txt
```

第 3/3 步：

```
neoHGT search -i <input_file.faa>
-m precomp -s <output_directory>/result.txt
-t <taxdump_dir>
-o <output_directory>
```

查看你的 ```<output_directory>``` 中的 ```.tsv``` 文件，如果它非空，就可以继续进行 [分析](#分析) 步骤。

___

如果由于某些原因 ```diamond``` 无法工作，备选方案是 ```blastp```，可在[这里](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)找到。但首先，你需要按照这篇[文章](https://www.ncbi.nlm.nih.gov/books/NBK569841/)来编译 blastp 数据库。

构建 blast 数据库的命令大致如下：

```
makeblastdb -parse_seqids -in db.faa
-out <blast_db_path>/<name_of_db>
-title <name_of_db>
-taxid_map <taxon.map>
```

这里的 taxon.map 文件可以通过运行以下命令获得：

```
gunzip -k taxon.map.gz
```

然后，输入以下命令：

```
blastp -query <input_file.faa>
-db <blast_db_path>/<name_of_db>
-outfmt "6 qaccver saccver pident evalue bitscore qcovhsp staxids"
> <output_directory>/result_blast.txt
```

最后，

```
neoHGT search -i <input_file.faa>
-m precomp -s <output_directory>/result_blast.txt
-t <taxdump_dir>
-o <output_directory>
```

查看你的 ```<output_directory>``` 中的 ```.tsv``` 文件，如果它非空，就可以继续进行 [分析](#分析) 步骤。

___

除了手动方式之外，下面这条命令会把搜索过程一站式完成，但不推荐使用，因为由于 ```diamond``` 的依赖版本升级，一个已知的兼容性 bug 已被发现并得到证实：

```
neoHGT search -i <input_file.faa.gz> -o . -m diamond -p 32 -d <diamond_db> -t <taxdump_dir>
```

要了解更多细节，请参阅 [HGTector 的 issue 页面](https://github.com/qiyunlab/HGTector/issues/114)。

## 分析

遵循上游 neoHGT，图形的输出格式是 `.pdf` 而非 `.png`，以获得更高的分辨率并便于发表。

分析所得 ```.tsv``` 文件的命令是：

```
neoHGT analyze -i <output.tsv>
-t <taxdump_dir>
-o <output_dir>
```

它会生成若干文件，其中最重要的是 ```scatter.pdf```，它是整条流水线的主要输出。

___

出于研究发表的严谨性考虑，你可能希望完全掌控诸如 ```--self-tax``` 和 ```--close-tax``` 这样的变量，否则程序会自动推断这些值。

```
neoHGT analyze -i <output.tsv>
-t <taxdump_dir>
-o <output_dir>
--self-tax <NCBI Taxonomy ID>
--close-tax <NCBI Taxonomy ID>
```
___

上面为了清晰起见，命令被拆分成多行，但在实际操作中你应该保持为一行，不要按 ```return``` 键。


## 卸载

要卸载，你可以运行：
```
conda remove -n neoHGT --all
```

注意：`bioconda` 是一个 conda *通道*，而不是一个包，因此没有什么可在那里卸载的；如果你还想移除 DIAMOND，可以像上面那样删除整个环境，或者运行 `conda remove -n neoHGT diamond`。


## 发表文献

上游项目 neoHGT（v1）是作者在 LOB - Laboratoire d'Optique et Biosciences 实习期间的成果，该团队在 [Nature Communications](https://www.nature.com/articles/s41467-023-36487-z) 上发表了一篇论文；neoHGT2GTDB 将这项工作作为其科学基线和引用来源。

> Filée, J., Becker, H. F., Mellottee, L., Eddine, R. Z., Li, Z., Yin, W., Lambry, J.-C., Liebl, U., & Myllykallio, H. (2023, February 15). Bacterial origins of thymidylate metabolism in Asgard Archaea and Eukarya. Nature Communications. Retrieved February 18, 2023, from https://www.nature.com/articles/s41467-023-36487-z



## 文献引用

> Zhu Q, Kosoy M, Dittmar K. HGTector: [an automated method facilitating genome-wide discovery of putative horizontal gene transfers](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-717). *BMC Genomics*. 2014. 15:717.


## 许可证

本仓库——非官方的 neoHGT2GTDB 适配——及其各上游项目均以 [BSD 3-clause](http://opensource.org/licenses/BSD-3-Clause) 许可证分发；参见完整的[声明](LICENSE)，其中载有每一条上游版权声明。

Copyright (c) 2026-, 曾子超 (Zeng Zichao, [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)) 以及 neoHGT2GTDB 的各位贡献者。

Copyright (c) 2023-, Zhihui LI.

Copyright (c) 2013-, [Qiyun Zhu](mailto:qiyunzhu@gmail.com) 和 [Katharina Dittmar](mailto:katharinad@gmail.com)——上游 HGTector 和 neoHGT 的作者，本项目正是从这两个项目衍生而来。
