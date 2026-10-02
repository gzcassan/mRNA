# B16-F10：从第一条 VCF 突变到候选筛选的实际证据

日期：2026-10-02；版本：0.1.0。

这次直接读取 OpenVax/Vaxrank 的 B16-F10 小鼠黑色素瘤测试文件，核对 BAM 与 SAM，检查 RNA reads，并调用 IEDB 的 NetMHCpan 4.1 BA 服务获得真实模型输出。这里的“真实”指实际文件和实际预测，不意味着已完成疫苗实验或临床验证。

**第一条突变的结论：Aldh1b1 的 G→C 在这个 RNA 子集中没有突变碱基支持，因此在本次采用的“必须有突变 RNA 支持”路线中被排除，不进入 mRNA 候选。** 为继续理解肽段与 MHC 分数，下文再追踪第五条 Wdr13 记录。

## 1. 我们实际拿到了什么

固定上游代码提交：`e108625abc3ec3e6408f9cb44c8be12710bb457b`。

- [b16.vcf：五条变异记录](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/data/b16.f10/b16.vcf)
- [b16.combined.bam：RNA 比对测试数据](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/data/b16.f10/b16.combined.bam)
- [b16.combined.sam：可直接阅读的文本形式](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/data/b16.f10/b16.combined.sam)

BAM 和 SAM 都包含 300 条比对记录。本次独立解码 BAM，并核对所有记录的名称、flag、染色体、起点、MAPQ、CIGAR 和碱基序列，确认与 SAM 一致；不是仅凭 SAM 文件名假定它们配套。

这是精选的测试子集，不是整份肿瘤 RNA 测序。文件校验值与上游 manifest 一致，但上游说明这些 B16 文件没有完整记录原始数据来源，因此不能宣称从原始研究 FASTQ 完成了全流程复现。[上游 TEST_DATA.md](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/TEST_DATA.md)

| 输入 | SHA-256 |
| --- | --- |
| b16.vcf | `d3d445b4e15f53c0861046f2ed819aa5cecab27e2401d5f8e1d734942513654a` |
| b16.combined.bam | `5c6161eca2410e8d8d54feee91ed3952ea04822face1d58be41593aaae0be21d` |
| b16.combined.sam | `b1dcf01943746f3d4861aa8d6f0ddcfbc7b39eaae9b57e93ae7d44aff40ceb37` |

上游演示脚本使用 `--mhc-predictor random`。随机分数能检查软件输出，不能回答哪些肽更可能结合 MHC。本次实际提交公开示例序列到 IEDB，使用固定方法 `netmhcpan_ba-4.1`。[上游演示脚本](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/run-vaxrank-b16-test-data.sh)、[IEDB API 文档](https://tools.iedb.org/main/tools-api/)

## 2. 第一条记录：逐列读 VCF

```text
##reference=mm10
#chr    pos         id  ref  alt  qual  filter  info
chr4    45802539    .   G    C    .     .       .
```

| 列 | 本条值 | 含义 |
| --- | --- | --- |
| CHROM | chr4 | 小鼠第 4 号染色体 |
| POS | 45802539 | mm10 参考体系中的位置，VCF 使用从 1 开始的坐标 |
| ID | . | 没有提供变异标识 |
| REF | G | 参考基因组此处是 G |
| ALT | C | 变异清单记录的替代碱基是 C |
| QUAL | . | 没有提供变异质量分数 |
| FILTER | . | 没有提供过滤状态；不能解释成明确的 PASS |
| INFO | . | 没有提供附加注释、深度或等位比例 |

这是简化测试 VCF：没有 FORMAT 和样本基因型列，也没有配对正常 DNA 证据。我们可以分析仓库提供的候选，但不能仅凭这行重新证明它是可靠的体细胞突变。

UCSC 的 mm10 `ncbiRefSeq` 注释将该位点定位到正链 Aldh1b1，转录本 `NM_028270.4`。基因名称是通过注释得到的，VCF 这行自身没有写基因名。

## 3. 找 RNA reads：有覆盖，为什么仍被淘汰？

我们不是看“基因表达量高不高”，而是看 **RNA 在这个具体位点有没有 C**。

一条实际记录：

```text
read:   HWI-D00273:119:C7FUMANXX:2:1302:1900:88137
FLAG:   147
RNAME:  chr4
POS:    45799162
MAPQ:   255
CIGAR:  8M3285N93M
局部:   TCTACAGCAGCT[G]CTCTCCCG
```

`8M3285N93M` 表示先比对 8 个碱基，跳过 3285 个参考碱基，再比对 93 个碱基。`N` 是 RNA 比对中可见的参考跳过区；必须按 CIGAR 映射，不能用“目标位置减起点”直接索引 read。

本条 read 的第二个比对块从 `45802455` 开始。目标 `45802539` 位于 read 的从 0 开始索引 `8 + (45802539−45802455) = 92`，也就是第 93 个碱基。读到的是 **G**，不是 C；该碱基的 Phred 质量为 35。MAPQ=255 表示此字段没有可用的常规质量值，不能把 255 当成最高置信度。

对整个测试文件进行同样检查：

| 计数单位 | 覆盖该碱基 | 支持 G | 支持 C |
| --- | ---: | ---: | ---: |
| 比对 reads | 14 | 14 | 0 |
| 按 read 名称合并的配对片段 | 13 | 13 | 0 |

配对片段计数避免把同一模板的两个 mates 重复计入；不是基于 UMI 的分子计数，也没有证明所有片段都来自独立原始分子。

这些结果与上游 [b16.not-expressed.vcf](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/data/b16.f10/b16.not-expressed.vcf) 的分类一致。但更精确的表述是“这个子集中没有突变 RNA 支持”，不是“Aldh1b1 基因完全不表达”。

```text
VCF 中有 G→C
    → RNA 确实覆盖该位点
    → 只看到 G，未看到 C
    → 无法建立 RNA 支持的突变蛋白片段
    → 本次 RNA 支持路线排除
    → 不进入 mRNA 候选
```

**这条记录没有继续计算突变肽的 MHC 分数。应写“不适用”，而不是“分数为 0”。** 即使可以从参考注释生成 DNA 假设肽，也不会补上缺失的突变 RNA 证据。其他允许 DNA-only 候选的配置应另行标记证据来源，不能冒充本路线。

## 4. 为继续看肽段：第五条 Wdr13 记录

```text
chrX    8125624    .    C    A    .    .    .
```

这次采用第五条，是因为第一条已经在 RNA 层被排除。第二条 Phip 也有 RNA 支持，但仓库提示其附近涉及另一处变化，需要处理同一密码子中的联合效应；不应把不同记录混成一条路线。

| Wdr13 RNA 证据 | 覆盖 | 参考 C | 替代 A |
| --- | ---: | ---: | ---: |
| 比对 reads | 73 | 48 | 25 |
| 按 read 名称合并的片段 | 65 | 44 | 21 |

例如实际 read `HWI-D00273:119:C7FUMANXX:2:2203:20863:75018`，起点 `8125526`，CIGAR 为 `101M`。目标位于索引 `8125624−8125526=98`：

```text
TCACAGTTGACG[A]TG
```

它支持 A，位点 Phred 质量为 38。以上 reads/片段计数也与上游 Wdr13 测试的预期一致；本次通过原始比对独立重算，没有将测试中的随机 MHC 分数拿来使用。[上游相关测试](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/test_mutant_protein_sequence.py)

## 5. 从 C→A 到 S460I：方向和阅读框都要对

UCSC mm10 注释中的 `NM_026137.5` 位于负链。我们按编码外显子拼接参考序列、取反向互补，再定位编码碱基，得到：

```text
基因组正向：C → A
编码方向：  G → T
编码位置：  c.1379G>T
密码子：    AGC → ATC
氨基酸：    S（丝氨酸）→ I（异亮氨酸）
蛋白位置：  p.S460I
```

这是按该 RefSeq 转录本得到的定位；其他转录本的蛋白编号应分别核对。RNA reads 提供序列证据，翻译依赖正确的参考阅读框，不是测到了真实蛋白。

围绕此位点取 25 aa 窗口，RNA 支持的序列为：

```text
KLQGHSAPVLDVIVNCDESLLASSD
            ^
            窗口第 13 位 I，对应蛋白第 460 位
```

所选编码窗口每个碱基至少由 13 条突变位点支持 reads 覆盖；这是逐碱基覆盖量，不等于 13 条 reads 各自覆盖整个窗口，也不是完整 Isovar 装配结果。

### 不能忽略的邻近变化

在支持目标 A 的 25 条 reads 中，全部也在 `chrX:8125622` 显示 C，mm10 参考为 A。这个邻近变化对应参考转录本的下一位 `F461V`。其体细胞或生殖系性质，不能靠这个 RNA 子集确定。

因此必须比较三个序列：

| 序列 | 25 aa 窗口 | 用途 |
| --- | --- | --- |
| RNA 观察背景 | `KLQGHSAPVLDVIVNCDESLLASSD` | 同时保留 I460 和 V461 |
| 仅恢复目标位点的比较序列 | `KLQGHSAPVLDVSVNCDESLLASSD` | 在 V461 背景中比较 S460 与 I460；是计算比较序列 |
| mm10 完整参考窗口 | `KLQGHSAPVLDVSFNCDESLLASSD` | 对照参考基因组，保留 S460 和 F461 |

第二行不能称为已经确认的正常样本序列。RNA 窗口另有远处参考差异，但不落在这个 25 aa 编码窗口内。

## 6. 切成肽段，并取得真实 MHC 模型输出

本次对上述三个 25 aa 序列扫描 8、9、10、11 aa 肽段，使用小鼠 `H-2-Kb` 和 `H-2-Db`。这不是人类 HLA 分型案例。

调用方法：IEDB API 的 `netmhcpan_ba-4.1`。完整请求保存在 `iedb-request.json`，原始响应保存在 `iedb-ba.tsv`。[IEDB API](https://tools.iedb.org/main/tools-api/)

每个窗口有 `18+17+16+15=66` 个连续子肽，两种 MHC 共 132 个配对。三个序列共返回 **396 条预测**。其中 RNA 观察序列包含目标 I460 的配对为 `(8+9+10+11)×2=76` 个。

| RNA 背景肽 | MHC | 预测 IC50，nM | BA percentile rank，% |
| --- | --- | ---: | ---: |
| SAPVLDVIV | H-2-Db | 797.39 | 0.20 |
| SAPVLDVI | H-2-Db | 2731.78 | 0.64 |
| HSAPVLDVI | H-2-Db | 5284.83 | 1.30 |
| SAPVLDVIV | H-2-Kb | 6236.89 | 6.20 |

`SAPVLDVIV` 是窗口第 6–14 位，I 位于肽的第 8 位。IC50 越低代表模型预测结合越强；rank 越低代表相对于模型背景肽更靠前。**0.20% 不是“0.20% 的免疫反应概率”。** 这是 BA 结合预测，不是 EL 呈递分数，也不是 T 细胞实验结果。[IEDB 结果解释](https://tools.iedb.org/mhci/help/)

比较同一个 9-mer 的三个背景：

| 序列来源 | 肽 | H-2-Db IC50，nM | BA rank，% |
| --- | --- | ---: | ---: |
| RNA 观察背景 | SAPVLDVIV | 797.39 | 0.20 |
| 仅恢复 I460 | SAPVLDVSV | 481.59 | 0.13 |
| mm10 参考 | SAPVLDVSF | 1226.08 | 0.29 |

因此不能笼统地说“突变比正常结合更强”：相对于只恢复 I460 的背景，结合变弱；相对于完整 mm10 参考，则变强。比较对象不同，结论就不同。结合差异仍不能证明 T 细胞特异性。

## 7. 它为什么进入或没进入候选？

必须先说明筛选规则。对 **包含目标 I460 的 76 个配对**：

| 明确规则 | 通过数量 | 判断 |
| --- | ---: | --- |
| 探索性 BA rank ≤2% | 3 | Wdr13 可继续作为探索候选 |
| IC50 <5000 nM | 2 | 有两项进入亲和力评分范围 |
| 严格 IC50 ≤500 nM | 0 | Wdr13 在这个目标位点规则下排除 |

≤2% 和 ≤500 nM 在这里是用于说明不同决策的显式规则，不是自动宣称的统一临床阈值。额外 F461V 来源的表位也没有在本表作为独立靶点评估。

为理解 Vaxrank 的一项经典评分，本次按上游测试参考公式手动重算单个窗口：

```text
expression_score = sqrt(支持突变的片段数) = sqrt(21) = 4.582576

单个配对评分（IC50 <5000 nM）：
score = (1 + exp(-350/150)) / (1 + exp((IC50-350)/150))

797.39 nM → 0.052893977
2731.78 nM → 0.000000139
target_epitope_score 合计 → 0.052894116

示例 combined_score = 4.582576 × 0.052894116 = 0.242391291
```

这是**固定 25 aa 窗口的手工公式重算**，不是完整 Vaxrank 排名，不是临床有效率。尚未执行全蛋白组自身序列过滤、全部窗口竞争、其他变异比较、连接区审查与构建容量选择。[上游评分参考](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/tests/_legacy_score_reference.py)、[表达评分实现](https://github.com/openvax/vaxrank/blob/e108625abc3ec3e6408f9cb44c8be12710bb457b/vaxrank/vaccine_peptide.py)

## 8. mRNA 终点：这次能下什么结论

| 对象 | 已完成的判断 | mRNA 状态 |
| --- | --- | --- |
| 第一条 Aldh1b1 G→C | 本子集 14 条覆盖 reads 全是 G；0 条支持 C | 本次 RNA 支持路线排除，不编码进候选 |
| 第五条 Wdr13 C→A | 有 RNA 支持；得到局部翻译序列和实际 MHC 预测 | 按宽松规则可继续；按 ≤500 nM 规则排除；最终构建未决定 |

这次没有生成或宣称已选出的 `full.fasta`。第一条的筛选过程已形成明确排除结果；第五条展示了实际序列如何产生分数，以及规则如何影响保留决定。完整 Vaxrank CLI 与 Isovar 装配尚未运行：本机 Linux 环境不可用，因此这次采用直接文件审计、RefSeq 阅读框核对和 IEDB 在线预测，不以随机模型或手工分数冒充全流程输出。

若采用探索性规则继续 Wdr13，需要先核实配对正常样本背景、检查全部参考蛋白中的自身匹配、对候选窗口及其他靶点统一排序，然后才决定是否加入 mRNA；即便加入计算构建，也仍需验证天然呈递、T 细胞识别与肿瘤功能。

## 9. 证据文件与复查方式

仓库中的 `results/b16-f10/` 保存：

- `Aldh1b1-reads.tsv`、`Wdr13-reads.tsv`：逐条 read 的起点、CIGAR、位点碱基、质量和局部序列。
- `summary.json`：原始计数、文件哈希与局部 RNA 窗口。
- `translation.json`：参考转录本、密码子、蛋白窗口与参考差异。
- `iedb-request.json`、`iedb-ba.tsv`：实际 API 输入和未改写响应。
- `Wdr13-peptide-comparisons.tsv`：76 个目标位点配对的三个背景比较。
- `ranking.json`：显式规则计数及单窗口评分。
- UCSC 查询响应与 `provenance.json`：参考序列、注释、查询参数和证据文件校验值。

两个复查脚本位于 `scripts/b16_trace/`。它们读取固定上游文件与保存的参考/API 响应，复算本笔记的结果，不运行完整 Vaxrank。Python 翻译使用 Biopython 1.88。

```powershell
git clone https://github.com/openvax/vaxrank.git work/vaxrank
git -C work/vaxrank checkout e108625abc3ec3e6408f9cb44c8be12710bb457b
python -m pip install biopython==1.88
python scripts/b16_trace/inspect_b16.py --data work/vaxrank/tests/data/b16.f10 --evidence results/b16-f10
python scripts/b16_trace/complete_b16_evidence.py --data work/vaxrank/tests/data/b16.f10 --evidence results/b16-f10
```

重新调用 IEDB 时，应使用保存请求中的方法、序列、MHC 和肽长；在线服务可能更新，结果差异需要另记版本与日期。保存的预测响应可供离线审计，不需要重复发送。

下一次完整运行的目标是：固定可工作的 Linux 环境、Isovar/Vaxrank/参考版本和真实预测器，完成所有候选的统一筛选，保留每条排除原因，再检查最终 manifest 中是否包含 Wdr13。
