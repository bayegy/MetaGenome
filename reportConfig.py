sum_steps_reads = """
			<p >1) 数据质控和去宿主序列：采用KneadData软件对原始数据进行质控(基于Trimmomatic)和去宿主(基于Bowtie2)，KneadData前和KneadData后，会用FastQC来检测质控合理性和效果。<sup>[7，8]</sup><br></p>
			<p >2）物种注释：使用Kraken2和自建的微生物数据库(从<a href="http://ccb.jhu.edu/software/kraken2/index.shtml?t=manual">Kraken官网</a>下载各个微生物数据库(细菌，真菌，古菌，病毒)，合并后再加入哥伦比亚大学实验室研究中新发现的一些细菌基因组数据）来鉴别样本中所含有的物种，再用Bracken来对样本中物种的实际相对丰度进行预测。相较于基于组装的物种注释，基于序列的宏基因组物种注释方法更加全面和准确。<sup>[9-12]</sup><br></p>
			<p >3) 常用功能数据库注释：从质控以及去除宿主基因的reads出发，使用HUMAnN2软件（基于DIAMOND），将各个样本的reads比对到数据库（UniRef90），根据UniRef90 ID 和各个数据库的对应关系，得到各个功能数据库的注释信息和相对丰度表。<sup>[19-22]</sup><br></p>
			<p >4) 基于物种丰度表和功能丰度表，可以进行丰度聚类分析，PCoA和NMDS降维分析（仅物种），样品聚类分析；当有分组信息时，可以进行LEfSe biomarker挖掘分析以及代谢通路比较分析，挖掘样品之间的物种组成和功能组成差异。<sup>[22，23]</sup><br></p>
			<p >5) 抗性基因注释：从去除宿主基因的clean reads出发，使用FMAP软件将各个样本的质控以及去除宿主基因的reads与抗生素抗性基因数据库CARD进行比对注释，可以获得抗性基因丰度分布情况。<sup>[21,24-25]</sup><br></p>
			<p >6) 另外，还可以基于标准分析结果，进行一系列高级信息分析（如 肠型分析，分箱分析，病原与宿主互作数据库(PHI)注释，分泌蛋白预测，III型分泌系统效应蛋白预测，细菌致病菌毒力因子(VFDB)注释，转移元件分析（MGE）等，更多详细信息请查看微生太宏基因组高级信息分析说明）；同时，结合环境因子、病理指标或特殊表型进行深入关联研究，能够为进一步深入研究和利用样品的物种和功能提供理论依据。<sup>[27-29]</sup>
			</p>
"""

sum_steps_asem = """
			<p >1) 数据质控和去宿主序列：采用KneadData软件对原始数据进行质控(基于Trimmomatic)和去宿主(基于Bowtie2)，KneadData前和KneadData后，会用FastQC来检测质控合理性和效果。<sup>[7，8]</sup><br></p>
			<p >2）物种注释：使用Kraken2和自建的微生物数据库(从<a href="http://ccb.jhu.edu/software/kraken2/index.shtml?t=manual">Kraken官网</a>下载各个微生物数据库(细菌，真菌，古菌，病毒)，合并后再加入哥伦比亚大学实验室研究中新发现的一些细菌基因组数据）来鉴别样本中所含有的物种，再用Bracken来对样本中物种的实际相对丰度进行预测。相较于基于组装的物种注释，基于序列的宏基因组物种注释方法更加全面和准确。<sup>[9-12]</sup><br></p>
			<p >3) 功能注释：运用MEGAHIT<sup>[46]</sup>软件，将所有样本去宿主基因后的clean reads进行组装，得到contigs; 运用prodigal软件，预测contigs中的基因序列; 再用Cd-hit软件，对得到的基因进行去冗余，得到去冗余基因; 使用Salmon软件，对去冗余基因进行定量; 使用eggnog-mapper, DIAMOND软件，对去冗余基因进行各个数据库的注释。统计各个数据库的基因相对丰度表。<sup>[19-22]</sup><br></p>
			<p >4) 抗性基因注释：运用DIAMOND软件，将去冗余基因比对到CARD数据库，得到CARD数据库的抗性基因注释信息，根据去冗余基因的丰度信息，统计抗性基因的相对丰度表。<sup>[24-25]</sup><br></p>
			<p >5) 基于物种丰度表和功能丰度表，可以进行丰度聚类分析，PCoA和NMDS降维分析（仅物种），样品聚类分析；当有分组信息时，可以进行LEfSe biomarker挖掘分析以及代谢通路比较分析，挖掘样品之间的物种组成和功能组成差异。<sup>[22，23]</sup><br></p>
			<p >6) 另外，还可以基于标准分析结果，进行一系列高级信息分析（如 肠型分析，分箱，病原与宿主互作数据库(PHI)注释，分泌蛋白预测，III型分泌系统效应蛋白预测，细菌致病菌毒力因子(VFDB)注释，转移元件分析（MGE）等，更多详细信息请查看微生太宏基因组高级信息分析说明）；同时，结合环境因子、病理指标或特殊表型进行深入关联研究，能够为进一步深入研究和利用样品的物种和功能提供理论依据。<sup>[27-29]</sup>
			</p>
"""

nav3_asem = """
					<a class="dropdown-item" href="#a3.2">3.2 组装</a>
					<a class="dropdown-item" href="#a3.3">3.3 基因预测，基因去冗余，基因定量</a>
"""

nav5_reads = """
					<a class="dropdown-item" href="#a5.1">5.1 KEGG数据库</a>
					<a class="dropdown-item" href="#a5.2">5.2 MetaCyc数据库</a>
					<a class="dropdown-item" href="#a5.3">5.3 EggNOG数据库</a>
					<a class="dropdown-item" href="#a5.4">5.4 GO数据库</a>
					<a class="dropdown-item" href="#a5.5">5.5 四级EC酶</a>
					<a class="dropdown-item" href="#a5.6">5.6 CAZy碳水化合物活性酶库</a>
"""


nav5_asem = """
					<a class="dropdown-item" href="#a5.1">5.1 KEGG数据库</a>
					<a class="dropdown-item" href="#a5.3">5.2 EggNOG数据库</a>
					<a class="dropdown-item" href="#a5.4">5.3 GO数据库</a>
					<a class="dropdown-item" href="#a5.6">5.4 CAZy碳水化合物活性酶库</a>
"""

find_gene_asem = """
			<div id="a3.2" class="anchor"></div>
			<h3>
				3.2 组装(<a href="../01-Assembly">../01-Assembly目录</a>)
			</h3>
			<p>组装步骤如下：</p>
			<ol>
				<li>经过质控和去宿主得到 Clean Data，使用MEGAHIT 组装软件进行组装分析( Assembly Analysis )；</li>
				<li>对于单个样品，组装时选取 MEGAHIT<sup>[46]</sup> 的默认参数进行组装，得到该样品的contigs；组装参数：--k-list 21,29,39,59,79,99,119,141 --min-contig-len 500</li>
				<li>将各样品质控和去宿主后的 Clean Data 采用 Bowtie2 软件比对至各样品组装后的contigs上，获取未被利用上的 PE reads；比对参数：--end-to-end, --sensitive</li>
				<li>将各样品未被利用上的 reads 放在一起，选取MEGAHIT 适合组装复杂的大型的宏基因组的参数进行组装，得到所有样品混合组装的contigs；组装参数：--presets meta-large --k-min 27 --k-max 127 --k-step 10 --min-contig-len 500</li>
				<li>对于单样品和混合组装生成的 contigs，过滤掉 500bp以下的片段，并进行统计分析和后续基因预测；</li>

			</ol>
			<div id="a3.3" class="anchor"></div>
			<h3>
				3.3 基因预测，基因去冗余，基因定量
			</h3>

			<p>运用prodigal<sup>[18]</sup>软件，预测组装得到的所有contigs中的基因序列；预测参数：-p meta（宏基因组模式）。</p>
			<p>使用Cd-hit<sup>[43]</sup>的默认参数，对prodigal预测得到的基因进行去冗余，得到去冗余基因；去冗余参数： -G 1（使用全局序列identity阈值） -c 0.9（默认的全局identity阈值）。</p>
			<p>使用Salmon<sup>[44]</sup>软件，将质控和去宿主后的Clean Data比对到去冗余基因上，从而计算去冗余基因的相对丰度RPM（reads per million）。Salmon定量参数: --meta (宏基因组模式)。</p>
			<p>使用emboss软件的transeq命令，将去冗余后的基因翻译为蛋白质序列，用于后续的比对和注释。</p>
"""

func_steps_asem = """
			<h3>基于组装的功能分析的基本步骤：</h3>

			<p>
				1）运用eggnog-mapper<sup>[45]</sup>软件（基于DIAMOND），将去冗余蛋白序列（对应去冗余基因核酸序列）比对到EggNOG数据库，得到蛋白的KEGG，GO，COG注释信息。比对参数：--seed_ortholog_evalue 0.00001<br>
			</p>
			<p>
				2）使用DIAMOND<sup>[26]</sup>，将去冗余蛋白序列比对到CAZy数据库，得到CAZy的注释信息。比对参数：-e 0.00001 （evalue 阈值） --id 80 （identity 阈值） --top 3（bit score 不低于最高分的3%）<br>
			</p>
			<p>
				3）根据去冗余基因的丰度表和各个数据库的注释信息，对于每个数据库，将注释到数据库相同基因家族(gene family)的去冗余基因丰度加和，筛除比对失败的去冗余基因，得到每个数据库基因家族的相对丰度表。<br>
			</p>

			<p>
				4）从各个数据库功能的相对丰度表出发，进行相对丰度柱形图展示，Circos图展示，丰度聚类热图展示，组间功能差异LEfSe分析，组间功能差异pair-wise多重比较DunnTest分析，KEGG通路图填色，功能与环境因子（或者其它组学数据）的相关性分析。<br>
			</p>
"""

func_steps_reads = """
			<h3>基于reads的功能分析基本步骤:</h3>
			<p>
				1）使用HUMAnN2软件<sup>[22]</sup>（2018年发表在Nature methods），将质控和去宿主之后的序列与蛋白质数据库（UniRef90）进行比对（基于DIAMOND）；<br>
			</p>
			<p>
				2）过滤掉比对失败的reads（HUMAnN2默认比对参数：translated_query_coverage_threshold = 90.0, prescreen_threshold = 0.01, evalue_threshold = 1.0, translated_subject_coverage_threshold = 50.0）；<br>
			</p>
			<p>
				3）统计UniRef90各个蛋白的相对丰度（RPKM ，reads per kilobase per million，校正样本比对成功reads(mapped reads)数以及基因长度后的丰度)。<br>
			</p>

			<p>
				4）根据UniRef90 的ID 和各个功能数据库ID的对应关系（主要来自<a href="https://www.genome.jp/linkdb/">LinkDB</a>），统计各个功能数据库对应功能相对丰度。<br>
			</p>
			<p>
				5）从各个数据库功能的相对丰度表出发，进行相对丰度柱形图展示，Circos图展示，丰度聚类热图展示，组间功能差异LEfSe分析，组间功能差异pair-wise多重比较DunnTest分析，显著差异功能物种来源柱形图分析，KEGG通路图填色，功能与环境因子（或者其它组学数据）的相关性分析。<br>
			</p>
"""

func_source_reads = """
			<p>
				在用LEfSe进行分组间显著性差异分析之后，我们也会将各个分组的特征功能挑出来，根据HUMAnN2的分析结果，得到功能的物种来源，绘制功能物种来源组成柱形图，如图5-4：
				<br>
			</p>
			<p align="center"><img src="FiguresTablesForReport/Figure5-4.png" alt="LEfSe未检测到特征KO"></p>
			<p align="center"><strong>图5-4 功能物种来源组成柱形图</strong></p>
			<dd align="center" style="font-size:80%;">说明：横坐标对应样本，以及样本分组，不同分组用不同颜色标出。纵坐标对应各个样本该功能的相对丰度，不同物种来源用不同颜色标出。</dd>
			<p>利用R语言的dunn.test包，我们对每一个功能的进行分组间两两多重比较，得到分析的表格，可在结果目录<a href="2-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/DunnTest/">2-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/DunnTest</a>查看分析结果。该分析首先用Kruskal-Wallis算出总的p值，再用Dunn.test核心算法进行多重比较，并用Bofferoni校正错误发现率。</p>

			<p>更多更详细的结果请查看结果目录<a href="2-FuctionAnalysis/1-KEGG/">2-FuctionAnalysis/1-KEGG</a><br></p>

			<div id="a5.2" class="anchor"></div>
			<h3 >5.2 MetaCyc数据库(<a href="2-FuctionAnalysis/2-Metacyc/">2-FuctionAnalysis/2-Metacyc</a>目录)</h3>
			<p ><a href="https://metacyc.org/">MetaCyc数据库</a>是一个阐明通过实验手段阐释代谢通路的数据库。MetaCyc目标是收集所有已知生命的代谢通路，是一个庞大而全面的数据库，目前包含了来自3009个不同生物的2722个代谢通路。MetaCyc的代谢网络包含了初生代谢，次生代谢，还包括相关的化合物、酶和基因。</p>

			<p >
				根据MetaCyc数据库的注释结果，绘制了各个通路在各样品的相对丰度统计图（图5-5, 相对丰度前20）。<br>
			</p>
			<p align="center"><img src="FiguresTablesForReport/Figure5-5.png" width="700" height="500"></p>
			<p align="center"><strong>图5-5 各个样品MetaCyc代谢通路相对分布情况柱形图（相对丰度前20的基因）</strong></p>
			<dd align="center" style="font-size:80%;">说明：图例中最多显示最优势的20个通路,余下的相对丰度较低的通路被归类为Other在图中展示。</dd>

			<p>更多更详细的结果请查看结果目录<a href="2-FuctionAnalysis/2-Metacyc/">2-FuctionAnalysis/2-Metacyc</a><br></p>

			<div class="alert alert-primary alert-dismissible">
				<button type="button" class="close" data-dismiss="alert">&times;</button>
				<strong>注意!</strong> 本结果中，所有的功能数据库注释结果的可视化方法基本一致，请参照其他功能数据库，抗性基因以及相关分析部分的可视化结果的解释，对该部分的结果进行解读。
			</div>
"""


func_source_asem = """
			<p>利用R语言的dunn.test包，我们对每一个功能的进行分组间两两多重比较，得到分析的表格，可在结果目录<a href="../../../2-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/DunnTest/">2-FuctionAnalysis/1-KEGG/4-SignificanceAnalysis/DunnTest</a>查看分析结果。该分析首先用Kruskal-Wallis算出总的p值，再用Dunn.test核心算法进行多重比较，并用Bofferoni校正错误发现率。</p>

			<p>更多更详细的结果请查看结果目录<a href="../../../2-FuctionAnalysis/1-KEGG/">2-FuctionAnalysis/1-KEGG</a><br></p>
"""


ec_doc_reads = """
			<div id="a5.5" class="anchor"></div>
			<h3>5.5 EC酶库(<a href="2-FuctionAnalysis/5-EC/">2-FuctionAnalysis/5-EC</a>)目录</h3>
			<p><a href="https://enzyme.expasy.org/">ENZYME</a>收录了酶的四级分类信息。EC编号或EC号是酶学委员会（Enzyme Commission）为酶所制作的一套编号分类法，每一个酶的编号都以字母“EC”起头，接着以四个号码来表示，这些号码代表逐步更细致的为酶作出分类。就如三肽胺基 蛋白酶的编号为EC3.4.11.4，当中的“EC3”是指水解酶（即以水来将分子分解的酶）；“EC3.4”是那些与肽键产生作用的水解酶；“EC3.4.11”是单指那些从多胜肽中分开胺基末端的水解酶；“EC3.4.11.4”则是从三肽中分开胺基末端的水解酶。</p>


			<p >
				在用LEfSe进行分组间显著性差异分析之后，我们也会将各个分组的特征功能挑出来，根据HUMAnN2的分析结果，得到功能的物种来源，绘制功能物种来源组成柱形图，如图5-8：
				<br>
			</p>

			<p align="center"><img src="FiguresTablesForReport/Figure5-8.png" alt="LEfSe未检测到特征酶活性基因"></p>
			<p align="center"><strong>图5-8 酶活性基因物种来源组成柱形图</strong></p>
			<dd align="center" style="font-size:80%;">说明：横坐标对应样本，以及样本分组，不同分组用不同颜色标出。纵坐标对应各个样本该功能的相对丰度，不同物种来源用不同颜色标出。</dd>

			<p>更多更详细的结果请查看结果目录<a href="2-FuctionAnalysis/5-EC/">2-FuctionAnalysis/5-EC</a><br></p>
			<div class="alert alert-primary alert-dismissible">
				<button type="button" class="close" data-dismiss="alert">&times;</button>
				<strong>注意!</strong> 本结果中，所有的功能数据库注释结果的可视化方法基本一致，请参照其他功能数据库，抗性基因以及相关分析部分的可视化结果的解释，对该部分的结果进行解读。
			</div>
"""
