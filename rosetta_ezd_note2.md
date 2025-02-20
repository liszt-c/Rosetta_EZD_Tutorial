#### 参考文献，2016 Rosetta and the Design of Ligand Binding Sites
#### 官方RosettaScripts教学: https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts  
#### 官方酶设计教程，主要看参数设置：https://docs.rosettacommons.org/docs/latest/application_documentation/design/enzyme-design  
#### 其他参考资料：https://zhuanlan.zhihu.com/p/70970820

# 1.1 relax protocol 松弛蛋白结构
>mpirun --oversubscribe -np 24 $ROSETTA3/bin/relax.static.linuxgccrelease -ignore_unrecognized_res -ignore_zero_occupancy false -use_input_sc -flip_HNQ -no_optH false -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -s MAL.pdb

ps: openmpi提示超过资源使用--oversubscribe来强制运行

生成的是MAL_0001.pdb，重命名为MAL_relax.pdb避免搞混了

# 1.2 对接，精确找到催化位置
用vina，柔性半柔性都可以，主要是为了定位。

# 1.3 生成200个小分子构象
注意仅支持sdf和mdl  
http://carbon.structbio.vanderbilt.edu/index.php/bclconf  
下载sdf文件名称改为SUR_conf.sdf
# 1.4 从小分子构象生成 params 文件
>$ROSETTA3/scripts/python/public/molfile_to_params.py --name=MOL --pdb=MOL --conformers-in-one-file MOL_conf.sdf 

其中，-n / --name 指定在 pdb 和 params 文件中用来表示配体名称的 3 字符缩写，需要注意的是，这里不能沿用 GLY 或者其他氨基酸、金属原子的缩写，否则生成的 params 文件会在后续突变结构时，与 Rosetta 自带的氨基酸或部分金属原子的 params 文件发生冲突；-p / --pdb 指定生成文件的命名。  
注意名称一致，不然后续会报错


# 2.1 定义约束文件(突变点位)
建立一个文本后缀.resfile，该文件定义参考https://zhuanlan.zhihu.com/p/70970820 ，
## 写入
>NATAA  
AUTO  
start  
1 X NATAA  
331 A NATAA  
73 A ALLAA  
360 A ALLAA  
361 A ALLAA  
389 A ALLAA  


注意，定多点突变而不是随机突变应该也可以写，看文章最后
### 参数
>ALLAA # 允许设计为20种氨基酸  
ALLAAxc # 允许设计为非半胱氨酸以外的所有氨基酸  
POLAR # 允许设计极性氨基酸(DEHKNQRST)  
APOLAR # 允许设计非极性氨基酸(ACFGILMPVWY)  
NOTAA <list of AAs> .. # 不允许设计为特定的氨基酸列表。(列表连续编写无空格)  
PIKAA <list of AAs> .. # 只允许设计为特定的氨基酸列表。(列表连续编写无空格)  
NATAA  # Repack当前氨基酸类型，只允许构象变化。  
NATRO  # NoRepack不允许构象变化。  
PROPERTY <property> # 只允许设计为有以下性质的氨基酸  

# 2.3 定义design.xml
详见文献

# 2.3运行RTOSETTA
### 参数
-ex1/-ex2 允许非天然氨基酸的侧链二面角扩展采样。  

-ex1aro/-ex2aro 专门针对芳香族残基（如Phe, Tyr）的扩展采样。  

level 4  最高采样级别，大幅增加尝试的旋转异构体数量，适用于高精度设计。  

-linmem_ig 10 
限制内部能量计算时使用的内存缓存大小。通过牺牲少量计算速度（约10%）来减少内存占用，避免大型任务（如多配体体系）崩溃。   

-restore_pre_talaris_2013_behavior
2013年Talaris力场更新后，部分能量计算方式改变。强制使用旧版能量函数逻辑，确保与历史脚本或特定研究的兼容性。

-parser:protocol design.xml 指定RosettaScripts XML流程文件  

-extra_res_fa MOL.params 添加非标准残基（如配体）的参数文件  

-s "MAL_relax.pdb MOL.pdb"  输入结构文件（蛋白+配体）  

-nstruct 5000 生成5000个独立结构

-out:file:scorefile design_results.sc 保存所有结构的评分结果
到design_results.sc


>mpirun --oversubscribe -np 12 $ROSETTA3/bin/rosetta_scripts.static.linuxgccrelease -ex1:level 4 -ex2:level 4 -ex1aro:level 4 -ex2aro:level 4 -linmem_ig 10 -restore_pre_talaris_2013_behavior -parser:protocol design.xml -extra_res_fa MOL.params -s "MAL_relax.pdb MOL.pdb" -nstruct 10000 -out:file:scorefile design_results.sc

ps: openmpi提示超过资源使用--oversubscribe来强制运行
### 48G内存开19个线程
>mpirun -np 19 $ROSETTA3/bin/rosetta_scripts.static.linuxgccrelease -ex1:level 4 -ex2:level 4 -ex1aro:level 4 -ex2aro:level 4 -linmem_ig 10 -restore_pre_talaris_2013_behavior -parser:protocol design.xml -extra_res_fa MOL.params -s "MAL_relax.pdb MOL.pdb" -nstruct 5000 -out:file:scorefile design_results.sc


# 3. 筛选
#### 参考文献，2016 Rosetta and the Design of Ligand Binding Sites
# 3.1 粗略筛选
###  3.1.1 准备一个文件(“metric_thresholds.txt”)，指定用于过滤设计运行输出的阈值。  
重要提示:阈值的精确值需要针对您的特定系统进行调整
>req total_score value < -1010 # measure of protein 
stability  
req if_X_fa_rep value < 1.0# measure of ligand 
clashes  
req ligand_is_touching_X value > 0.5# 1.0 if ligand 
is in pocket  
output sortmin interface_delta_X# binding energy  

### 3.1.2 使用rosetta自带脚本DesignSelect.pl筛选对接运行的初始指标。
将生成一个文件(“filtered_pdbs.txt”)，其中包含通过指标截止值的输出pdb列表。
PS: 这个命令在使用时删除design_results.sc前面的SEQUENCE，不然会报错 
>perl ${ROSETTA3}/src/apps/public/enzdes/DesignSelect.pl -d design_results.sc -c metric_thresholds.txt -tag_column last > filtered_designs.sc

-tag_column last: 标签（如PDB文件名）位于每行的最后一列

>awk '{print $NF ".pdb"}' filtered_designs.sc > filtered_pdbs.txt

# 3.2 附加指标的计算与筛选
### 3.2.1 重新计算附加指标
最好只在一组已经预先筛选过的结构上运行。生成指标后，可以按照前两个步骤中的方式过滤结构。将生成一个分数文件(“design_interfaces.sc”)，其中包含所选pdb的计算指标值。

>${ROSETTA}/main/source/bin/InterfaceAnalyzer.linuxgccrelease -interface A_X -compute_packstat -pack_separated -score:weights ligandprime -no_nstruct_label -out:file:score_only design_interfaces.sc -l filtered_pdbs.txt -extra_res_fa LIG.params
# 3.2.2筛选附加指标
这些命令类似于步骤3.1中使用的命令，但针对的是新的design_interfaces.sc分数文件，并且具有新的阈值文件。

阈值文件类格式如下，和前述一样，也需要根据体系调整(文件名可为metric_thresholds2.txt)：
>req packstat value > 0.55 # packing metric; 0-1 higher better  
req sc_value value > 0.45# shape complementarity; 0-1 higher better  
req delta_unsatHbonds value < 1.5# unsatisfied hydrogen bonds on binding  
req dG_separated/dSASAx100 value < -0.5 # binding energy per contact area  
output sortmin dG_separated# binding energy

### 调用 DesignSelect.pl
>perl ${ROSETTA}/main/source/src/apps/public/enzdes/DesignSelect.pl -d design_interfaces.sc -c metric_thresholds2.txt -tag_column last > filtered_design_interfaces.sc

>awk '{print $NF ".pdb"}' filtered_design_interfaces.sc > filtered_pdbs.txt

