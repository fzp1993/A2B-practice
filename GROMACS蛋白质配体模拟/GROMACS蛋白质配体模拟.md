# GROMACS蛋白质配体模拟

## 数据准备

我们关注的是7AIT和8UK的结合，因此我们需要从[RCSB](https://www.rcsb.org/ "RCSB")下载蛋白质pdb文件，这里选择PDB ID是7AIT的蛋白质。可以从download下载全部的蛋白质文件，也可以直接搜索PDB ID下载

<p align="center"><img src="pic/下载3HTB.png" alt="下载3HTB" width="50%"/></p>

7AIT.pdb中包含了3HTB和它所能结合的所有配体的信息，其中包括8UK的信息。由于GROMACS提供的所有力场都不能识别8UK配体，因此我们需要分别来准备蛋白质和配体的拓扑文件。

## 准备蛋白质拓扑文件

### 准备蛋白质pdb文件

* 在准备蛋白质文件之前，我们需要先将pdb文件中关于8UK的信息提取出来。进入ubuntu系统，输入命令保存8UK.pdb文件
  
  ```
  grep 8UK 7ait.pdb > 8uk.pdb
  ```

* 用PyMOL删除结晶水和配体信息。但是不一定所有的蛋白质处理过程都会删除结晶水，需要根据要求来进行处理。
  
  1、在命令框中输入“cd E:\pymol\7AIT_8UK”跳转到pdb文件所在文件夹。

  2、从File中选Open，打开pdb文件。

  3、红色的十字标就是结晶水。

  4、从A中选择remove waters来删除结晶水。

  5、从命令框中输入“indicate hetatm”，选出全部的hetatm，然后选择remove atoms，删除hetatm。其实结晶水也包括在hetatm中。

  6、从File中选择Expert Molecule来保存文件，保存类型要选择PDB。

  <p align="center"><img src="pic/原始7AIT.png" alt="原始7AIT" width="50%"/></p>
  <p align="center"><img src="pic/处理后的7AIT.png" alt="处理后的7AIT" width="50%"/></p>

* 有时候pdb文件会有氨基酸原子缺失，必须要将他们补充上。使用SPDBV打开pdb文件就可以直接补充缺失的原子。然后选择保存（File->Save->Current Layer）即可。

### 准备力场文件

* 在GROMACS中，已经封装了许多力场，根据需要去选择合适的力场。由于蛋白质配体结合更适合用CHARMM36力场，因此需要从[CHARMM36官网](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs "CHARMM36官网")下载最新的CHARMM36力场。

* 下载力场压缩包后，放到GROMACS的力场文件夹中，一般是在“/usr/local/gromacs/share/gromacs/top”里面（找不到也没关系，直接在下一步生成拓扑文件时，会让你选择力场文件，这时他会告诉你力场文件的文件夹位置。然后“ctrl+c”退出来，进入力场文件夹添加力场即可）。在ubuntu中执行命令安装
  
  ```
  tar -zxvf charmm36-jul2022.ff.tgz
  ```

### 生成蛋白质拓扑文件

* 输入命令生成蛋白质拓朴文件
  
  ```
  gmx pdb2gmx -f 7ait_clean.pdb -o 7ait_processed.gro
  ```

  需要注意的是，可能会显示“-bash: gmx: command not found”。这时候就需要激活一下GROMACS

  ```
  source /usr/local/gromacs/bin/GMXRC
  ```

* 选择力场
  
  ```
  Select the Force Field:
  From '/usr/local/gromacs/share/gromacs/top':
  1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
  2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
  3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
  4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
  5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
  6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
  7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
  8: CHARMM all-atom force field
  9: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
  10: GROMOS96 43a1 force field
  11: GROMOS96 43a2 force field (improved alkane dihedrals)
  12: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
  13: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
  14: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
  15: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
  16: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
  ```

  这里第8个就是我们添加的力场，输入8即可。

* 选择水模型
  
  这里选择默认的水模型即可，也就是CHARMM的TIP3P

  ```
  Select the Water Model:
  1: TIP3P      CHARMM-modified TIP3P water model (recommended over original TIP3P)
  2: TIP3P_ORIGINAL Original TIP3P water model
  3: SPC        SPC water model
  4: SPCE       SPC/E water model
  5: TIP5P      TIP5P water model
  6: TIP4P      TIP4P water model
  7: TIP4PEW    TIP4P/Ew water model
  8: None
  ```

  最后生成三个文件7ait_processed.gro、posre.itp、topol.top。

## 准备配体分子拓扑文件

实际上，可以从RCSB中下载到配体分子的mol2文件

### 添加氢原子

* 由于一般晶体结构中不包含氢原子，但是CHARMM是全原子力场，因此需要补充氢原子。使用[Avogadro软件](https://sourceforge.net/projects/avogadro/ "Avogadro软件")Avogadro软件来添加氢原子。

* 进入Avogadro，打开之前做好的8uk.pdb，然后选择“Buid->Add Hydrogens”。最后保存为mol2文件。
  
  <p align="center"><img src="pic/原始8UK.png" alt="原始8UK" width="50%"/></p>
  <p align="center"><img src="pic/加氢后的8UK.png" alt="加氢后的8UK" width="50%"/></p>

### 修改mol2文件

* 生成的8uk.mol2文件需要对几个地方进行修改，首先要用"8UK"代替原来的“*****”
  
  ```
  @<TRIPOS>MOLECULE
  *****
  ```

  ```
  @<TRIPOS>MOLECULE
  8UK
  ```

* 修改残基名称和残基序号，例如
  
  ```
  ...
  78 CAA        -9.9470  112.4280   77.8010 C.3   601  8UK      0.0790
  79 H          -0.3904   48.1268   53.2378 H       1  UNL1        0.0334
  ...
  ```

  ```
  ...
  78 CAA        -9.9470  112.4280   77.8010 C.3   601  8UK      0.0790
  79 H          -0.3904   48.1268   53.2378 H       1  8UK        0.0334
  ...
  ```

* 使用[sort_mol2_bonds.pl](http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.pl "sort_mol2_bonds.pl")对“@<TRIPOS>BOND”的键进行升序排列。在ubuntu中执行命令
  
  ```
  perl /mnt/GROMACS/sort_mol2_bonds.pl 8uk.mol2 8uk_fix.mol2
  ```

* 使用[CGenFF](https://cgenff.umaryland.edu/initguess/ "CGenFF")生成配体分子的拓扑文件，CGenFF是一个在线网站，直接上传mol2文件就可以生成CHARMM的"stream"文件。需要注意的是，这里使用的8uk_fix.mol2中的氢原子名称都是H，因此需要先把8uk_fix.mol2的氢原子H后面加上数字，然后再用CGenff生成拓扑文件。

* 使用cgenff_charmm2gmx_py3_nx2.py脚本将8uk.str转换为GROMACS可识别的文件，要注意的是需要先安装networkx 2.3。
  
  ```
  python /mnt/GROMACS/cgenff_charmm2gmx_py3_nx2.py 8UK 8uk_fix.mol2 8uk.str /usr/local/gromacs/share/gromacs/top/charmm36-jul2022.ff
  ```

  输出下面的信息后，说明转换成功

  ```
  ============ DONE ============
  Conversion complete.
  The molecule topology has been written to 8uk.itp
  Additional parameters needed by the molecule are written to 8uk.prm, which needs to be included in the system .top
  
  PLEASE NOTE: If your topology has lone pairs, you must use GROMACS version 2020 or newer to use 2fd construction
  Older GROMACS versions WILL NOT WORK as they do not support 2fd virtual site construction
  ============ DONE ============
  ```

  一共输出四个文件：8uk.top、8uk.itp、8uk_ini.pdb、8uk.prm。

## 构建蛋白质配体复合物

* 将8uk_ini.pdb转为gro格式，并和7ait_processed.gro进行拼接，最后修改原子个数为16851
  
  ```
  gmx editconf -f 8uk_ini.pdb -o 8uk.gro
  ```

  ```
  ...
  535THR    OT116694   0.052  12.224   5.762
  535THR    OT216695  -0.175  12.289   5.706
  18UK   CAU1    1  -0.056   4.817   5.429
  18UK   CBL1    2   0.007   4.945   5.484
  ...
  ```

* 构建拓扑文件只需要修改topol.top即可，在最后添加“#include "8uk.itp"”，在顶端添加8uk.prm的二面角信息，同时需要在[ molecules ]中加入配体分子。
  
  ```
  ; Include ligand topology
  #include "8uk.itp"
  
  ; Include water topology
  #include "charmm36-jul2022.ff/tip3p.itp"
  ```

  ```
  ; Include forcefield parameters
  #include "charmm36-jul2022.ff/forcefield.itp"
  
  ; Include ligand parameters
  #include "8uk.prm"
  ```

  ```
  [ molecules ]
  ; Compound        #mols
  Protein_chain_A     1
  Protein_chain_B     1
  8UK                 1
  ```

## 溶剂化

* 定义模拟盒子并添加水
  
  ```
  gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
  gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
  ```

## 添加离子

* 体系中溶剂化的蛋白质是带电荷的，从topol_Protein_chain_A.itp和topol_Protein_chain_B.itp文件中的[ atoms ]中的最后一行可以看到qtot -7和-8，说明体系中总电荷是-15，需要加入离子平衡电荷。首先下载[ions.mdp](http://www.mdtutorials.com/gmx/complex/Files/ions.mdp "ions.mdp")文件，然后执行命令
  
  ```
  gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
  ```

* 把得到的tpr文件传递给genion
  
  ```
  gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
  ```

  这时会问我们用什么分子来替换成钠离子，我们选择13，也就是水分子

  ```
  Will try to add 15 NA ions and 0 CL ions.
  Select a continuous group of solvent molecules
  Group     0 (         System) has 231792 elements
  Group     1 (        Protein) has 16695 elements
  Group     2 (      Protein-H) has  8490 elements
  Group     3 (        C-alpha) has  1064 elements
  Group     4 (       Backbone) has  3192 elements
  Group     5 (      MainChain) has  4254 elements
  Group     6 (   MainChain+Cb) has  5230 elements
  Group     7 (    MainChain+H) has  5258 elements
  Group     8 (      SideChain) has 11437 elements
  Group     9 (    SideChain-H) has  4236 elements
  Group    10 (    Prot-Masses) has 16695 elements
  Group    11 (    non-Protein) has 215097 elements
  Group    12 (          Other) has   156 elements
  Group    13 (            8UK) has   156 elements
  Group    14 (          Water) has 214941 elements
  Group    15 (            SOL) has 214941 elements
  Group    16 (      non-Water) has 16851 elements
  Select a group: 15
  ```

## 能量最小化

* 对体系solv_ions.gro进行能量最小化设置，首先下载[em.mdp](http://www.mdtutorials.com/gmx/complex/Files/em.mdp "em.mdp")文件，然后执行命令
  
  ```
  gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
  ```

  注意：有的时候添加了氯离子，会运行报错 “No such moleculetype CL”。解决方法是将topol.top文件和solv_ions.gro文件中[ molecules ]部分下的CL改为CLA。同时需要注意solv_ions.gro的缩进要保持格式规范。

* 然后调用mdrun进行能力最小化
  
  ```
  gmx mdrun -v -deffnm em
  ```

  如果出现下面的错误，就在命令后面加上 “-ntmpi 1”

  ```
  Fatal error:
  There is no domain decomposition for 8 ranks that is compatible with the given
  box and a minimum cell size of 6.11243 nm
  Change the number of ranks or mdrun option -rdd or -dds
  Look in the log file for details on the domain decomposition
  ```

## 平衡

### 约束配体

* 给配体创建一个包含氢原子之外的所有原子索引组
  
  ```
  gmx make_ndx -f 8uk.gro -o index_8uk.ndx
  ...
   > 0 & ! a H*
   > q
  ```

* 执行genrestr模块，然后选择刚才得到的索引组，是group 3
  
  ```
  gmx genrestr -f 8uk.gro -n index_8uk.ndx -o posre_8uk.itp -fc 1000 1000 1000
  ```

* 把信息添加到拓扑文件中
  
  ```
  ; Include Position restraint file
  #ifdef POSRES
  #include "posre.itp"
  #endif
  
  ; Include ligand topology
  #include "8uk.itp"
  
  ; Ligand position restraints
  #ifdef POSRES
  #include "posre_8uk.itp"
  #endif
  
  ; Include water topology
  #include "./charmm36-mar2019.ff/tip3p.itp"
  ```

### 热浴

* 控制温度耦合，但是不要单独耦合体系中的每一种分子。首先创建一个组合了蛋白质和配体的特殊索引组。可以通过make_ndx处理任意一个包含完整体系的坐标文件来实现，并合并蛋白质和8UK.
  
  ```
  gmx make_ndx -f em.gro -o index.ndx
  ...
   > 1 | 13
   > q
  ```

* 在[nvt.mdp](http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp "nvt.mdp")设置"tc-grps = Protein_8UK Water_and_ions"，并执行NVT平衡。
  
  ```
  gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
  gmx mdrun -deffnm nvt
  ```

* 用[npt.mdp](http://www.mdtutorials.com/gmx/complex/Files/npt.mdp "npt.mdp")执行NPT平衡。
  
  ```
  gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
  gmx mdrun -deffnm npt
  ```

## MD模拟

用[md.mdp](http://www.mdtutorials.com/gmx/complex/Files/md.mdp "md.mdp")执行NPT平衡。

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10
```

## PyMOL可视化轨迹

* 先将复合物放到之前定义的盒子中间，并选择Protein作为中心，并以Sytem作为输出
  
  ```
  gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
  ```

  执行旋转和平移拟合来获得更平滑的可视效果，并选择Backbone执行最小二乘法拟合蛋白质的骨架，然后选择System输出

  ```
  gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans
  ```
  
* 获取轨迹的第一帧
  
  ```
  gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o start.pdb -dump 0
  ```

* 打开PyMOL，加载start.pdb，然后加载轨迹文件md_0_10_center.xtc
  
  ```
  load_traj md_0_10_center.xtc,start=1,stop=10000000,interval=20
  ```

* 大功告成
  
  <p align="center"><img src="pic/complex.gif" alt="复合物轨迹动画" width="80%"/></p>