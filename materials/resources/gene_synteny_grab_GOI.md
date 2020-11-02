### Grab Gene of Interest

Below is the example documentation for how the file of genes that we were interested in (`bin2_cys.txt`) was generated. 

Load annotation files in R and find genes of interest. 

```R
library("tidyverse")
library("dplyr")

dataFiles <- lapply(Sys.glob("*.aa"), read.delim,sep = "\t")
head(dataFiles[[1]])
```

<table>
<caption>A data.frame: 6 × 33</caption>
<thead>
	<tr><th scope=col>Query.Gene</th><th scope=col>GC.</th><th scope=col>Contig.name</th><th scope=col>Start.position</th><th scope=col>Stop.position</th><th scope=col>Orientation</th><th scope=col>Query.sequence</th><th scope=col>Signalling</th><th scope=col>Signal.confidence</th><th scope=col>Target.gene..UniProt.</th><th scope=col>⋯</th><th scope=col>Coverage.2</th><th scope=col>E.value.2</th><th scope=col>Description.2</th><th scope=col>Taxonomy.2</th><th scope=col>Target.gene..Pfam.</th><th scope=col>E.value.3</th><th scope=col>Description.3</th><th scope=col>Target.gene..TIGRfam.</th><th scope=col>E.value.4</th><th scope=col>Description.4</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>⋯</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>caef037ad318582c885d129aee6008b1_3042</td><td>69.01</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MSPDSPHARLLRLATLASLGVAVTLVIGKAIAWWLSGSVSLLAGLTDSLLDSAASLLNLLAVHYSLRPADDDHRYGHGKAEALAGLAQAMFIGASAVLVAIQAVEQLRDPQPLGAATWGIAVMLLSLLLTALLLAFQTHVIRRTGSTAIRADSLHYRSDLLLNASILVALLLARYGWPQSDGLFGLGIALYILWSALQIARESVSILMDRELPTDISERMRELACSVPGVLGAHDVRTRISGNHWFVQLHLELPGELSLSVAHALTDRAVLAIKREYPKAEVLIHADPQEVVGKETVS*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>-</td><td>-</td><td>A0A069Q3E6</td><td>⋯</td><td>99.664430</td><td>-1.000000</td><td>transporter                                 </td><td>Pseudomonas aeruginosa PAO1</td><td>PF01545.21; PF16916.5                         </td><td>0.000000; 0.000000                    </td><td>Cation_efflux; ZT_dimer      </td><td>TIGR01297           </td><td>0.000000          </td><td>TIGR01297           </td></tr>
	<tr><td>caef037ad318582c885d129aee6008b1_3043</td><td>73.33</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MNGGWWRMLGVGALSLSLCVQAQADDGRGGRGGPGGERGVAHGQPGPGRGEAYGGSWGPHPDWRPGRVVEALPRGHLRVPYRGGEYFFHDGYWYRPDGPRYVLVVPPRGVRVRSLPPYAEQVWLGSVLYFLAAGTYYLWHADTREYEVVSPPPRAAGPGYPVSVPATGYDVVAYPARGQGPDQQSRDRYECHRWAVGESGFDPAGATYTPAAEVANGYRRAMAACLSGRGYSVN*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>-</td><td>-</td><td>A0A335Q0S2</td><td>⋯</td><td>99.559471</td><td>-1.000000</td><td>hypothetical protein                        </td><td>Pseudomonas aeruginosa PAO1</td><td>Not annotated                                 </td><td>-                                     </td><td>-                            </td><td>Not annotated       </td><td>-                 </td><td>-                   </td></tr>
	<tr><td>caef037ad318582c885d129aee6008b1_3040</td><td>71.28</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MISLPIDAVVPALRQALGAQHQAVLEAPPGAGKTTRVPLALLDEPWLAGQRILMLEPRRLAARAAAERLAAELGEKVGETVGYRIRLESRVGPKTRIEVVTEGILARRLQDDPALDGVGLVIFDEFHERSLDADLALALTLNGRELLRDEPPLKVLVMSATLEGERLAALLGEAPVVRSEGRMFPVDIRWGRPAQPGEFIEPRVQQAVLQALAEESGSVLVFLPGQAEIRRVHEGLREALGGRPEVLLCPLHGELDLAAQRAAIEPASRGTRKVVLATNIAETSLTIDGVRVVIDAGLARVPRFDPGSGMTRLETQRISRASATQRAGRAGRLEPGVCYRLWSESQHEQLPAYGTAEILQADLAGLALQLARWGVAPEELAWLDAPPAAAYAQARELLGRLGALNASGALSAHGQAMAELPTHPRIAHLLLRGQALGLGELACDVAALLGERDIQRGGGADLHSRLALLAGEARTGASRGAVQRARQLARQFRGYLRGAASEAVVDPGHPRWLGCLLAFAYPDRIARQRRAGGGDYRLANGRAAQFGEPDSLMKQPWLVIADLGSRQGQREERIYLAAELDPRLFDTVLAEQVSQRDELQWDEREGVLRAERQRRVGELVLSSEALPGLDEAARSQALLGLVRRKGLELLPWTPELRQWQARIGLLRRLDLEDKGESEWPDVSDAALLERLEEWLPAYLGKVTRLAHFANLDLASILAGLLPWPLPQRLDEWAPKTLEVPSGSRIRLDYSETPPILAVRLQELFGLGDTPRIAQGRLAVKLHLLSPAHRPVQVTQDLANFWRSTYAEVKKDLKGRYPKHYWPDDPLVAEATARAKPRK*</td><td>-</td><td>-</td><td>A0A219TCU6</td><td>⋯</td><td>99.880668</td><td>-1.000000</td><td>K03579: ATP-dependent helicase              </td><td>Pseudomonas aeruginosa PAO1</td><td>PF00270.29; PF04408.23; PF00271.31; PF08482.10</td><td>0.000000; 0.000000; 0.000000; 0.000000</td><td>DEAD; HA2; Helicase_C; HrpB_C</td><td>TIGR01967; TIGR01970</td><td>0.000000; 0.000000</td><td>TIGR01967; TIGR01970</td></tr>
	<tr><td>caef037ad318582c885d129aee6008b1_5549</td><td>67.10</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MRLDRFLANLPELSRRDAQMLIAGGRLRVDGQVVRDGTHEVRVFSRVELDERLLQAGKPARYYMLHKPMGCVSATRDPQHPTVLDLFPEDLRDDLHIGGRLDYNTTGLMLLTNDGQWSRRLTLPGSRRDKVYYVETEAPIDQRYVDAFAAGLHFRFEDLVTLPAQLEPIGPRCARLTLHEGRYHQVKRMFGHFDNKVLRLHRERMGDICLDPQLPPGAWRPLNAEEICSV*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>-</td><td>-</td><td>A0A2X2DA31</td><td>⋯</td><td>99.565217</td><td>-1.000000</td><td>K06183: 16S rRNA pseudouridine(516) synthase</td><td>Pseudomonas aeruginosa PAO1</td><td>PF00849.22; PF01479.25                        </td><td>0.000000; 0.000000                    </td><td>PseudoU_synth_2; S4          </td><td>TIGR00005; TIGR00093</td><td>0.000037; 0.000000</td><td>TIGR00005; TIGR00093</td></tr>
	<tr><td>caef037ad318582c885d129aee6008b1_5548</td><td>65.61</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MKRSLAVLLLAFAVGSLHAAPASPANAPLPPADIKARLKSMKPGEYLWYPEVSPQGPVTIVVSLTEQKAYIYRNGIAIGVSTLSSGKKGRETPTGVFSILQKSVDHKSDLYNSAPMPYMQRLTWDGIALHAGNLPGYPASHGCIRLPMAFAKKLYGITGFSSTTVIISNASSAPKEVDHPGLLAPTVADGRPVTLPVEPGEIAWNDPGAERGPLSVLISRADQRAYVYRGGQRIGTAPVQLPARPQTGMAVFSLLEKPTPEEISETSPNLRWSVVQVSNPGQGLTPSEQLGKIRMDPTFVREMLGAMDVGSTLVITDWASTRDTHSDSDFTVIATETPNNPTKSLKK*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           </td><td>-</td><td>-</td><td>A0A2C9WZZ2</td><td>⋯</td><td>99.711816</td><td>-1.000000</td><td>hypothetical protein                        </td><td>Pseudomonas aeruginosa PAO1</td><td>PF03734.14                                    </td><td>0.000000                              </td><td>YkuD                         </td><td>Not annotated       </td><td>-                 </td><td>-                   </td></tr>
	<tr><td>caef037ad318582c885d129aee6008b1_3041</td><td>66.67</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MRLATRVSCAVLSAALLSQLSACGTLFYPERRGQIDGRIDPAIVAFDAIGLLFYIIPGLIAFGVDFATGAIYLPDAKYSVAPEVLKDAVGADGKVDNRKLKAILEREIGRELPLDDPRLIRHSGSVEQLAAYGLKPAA*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>-</td><td>-</td><td>A0A072ZQ61</td><td>⋯</td><td>99.275362</td><td>-1.000000</td><td>hypothetical protein                        </td><td>Pseudomonas aeruginosa PAO1</td><td>Not annotated                                 </td><td>-                                     </td><td>-                            </td><td>Not annotated       </td><td>-                 </td><td>-                   </td></tr>
</tbody>
</table>

Here, we are looking for the gene cluster that consists of *sbp* (sulfur binding protein) and *cysUWA* using KEGG databases (thier KO numbers are: K02048,K02046,K02047,K02045)

```R
GOI <- data.frame()

for (keggid in c("K02048","K02046","K02047","K02045")) {
    GOI <- rbind(GOI,dataFiles[[1]][grep(keggid,dataFiles[[1]]$Description.2),])}

GOI
```

<table>
<caption>A data.frame: 5 × 33</caption>
<thead>
	<tr><th></th><th scope=col>Query.Gene</th><th scope=col>GC.</th><th scope=col>Contig.name</th><th scope=col>Start.position</th><th scope=col>Stop.position</th><th scope=col>Orientation</th><th scope=col>Query.sequence</th><th scope=col>Signalling</th><th scope=col>Signal.confidence</th><th scope=col>Target.gene..UniProt.</th><th scope=col>⋯</th><th scope=col>Coverage.2</th><th scope=col>E.value.2</th><th scope=col>Description.2</th><th scope=col>Taxonomy.2</th><th scope=col>Target.gene..Pfam.</th><th scope=col>E.value.3</th><th scope=col>Description.3</th><th scope=col>Target.gene..TIGRfam.</th><th scope=col>E.value.4</th><th scope=col>Description.4</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>⋯</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1286</th><td>caef037ad318582c885d129aee6008b1_465 </td><td>64.46</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MKRLFSASLLAAGLALGGAAHAAQPLLNVSYDVMRDFYKEYNPAFQKYWKAEKGENITIQMSHGGSSKQARSVIDGLPADVITMNQATDIDALADNGGLVPKDWATRLPNNSAPFTSATVFIVRKGNPKALKDWPDLLKDGVQVVVPNPKTSGNGRYTYLSAWGYVLKNGGDENKAKEFVGKLFKQVPVLDTGGRAATTTFMQNQIGDVLVTFENEAEMIAREFGRGGFEVVYPSVSAEAEPPVAVVDKVVEKKGSRAQAEAYLKYLWSDEGQTIAANNYLRPRNPEILAKFADRFPKVDFFSVEKTFGDWRSVQKTHFIDGGVFDQIYSSN*</td><td>-</td><td>-</td><td>A0A069Q3Z9</td><td>⋯</td><td>99.698795</td><td>-1.000000</td><td>K02048: cysP; sulfate ABC transporter substrate-binding protein           </td><td>Pseudomonas aeruginosa PAO1</td><td>PF13531.6; PF13343.6                                                                                              </td><td>0.000000; 0.000000                                                                                </td><td>SBP_bac_11; SBP_bac_6                                                           </td><td>TIGR00971; TIGR01256; TIGR03261                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   </td><td>0.000000; 0.000011; 0.000042                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          </td><td>TIGR00971; TIGR01256; TIGR03261                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   </td></tr>
	<tr><th scope=row>3620</th><td>caef037ad318582c885d129aee6008b1_3761</td><td>65.47</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MSLRRFALAALASALVSGPVAAATQLLNVSYDPTRELYQAYNAAFIKHWKAQGGEDLTVQQSHGGSGKQARAVIDGLKADVVTLALAGDIDELHKLGKLLPADWQARLPDNSTPYTSTIVFLVRKGNPKGIKDWGDLTKEGVEVITPNPKTSGGARWNFLAAWAWAKKQYGSDEKAKDYVQALYKHVPVLDTGARGSTITFVNNQIGDVLLAWENEAFLAKKEQGGENFEIVVPSISILAEPPVAVVDKVVEKKGTRKVAEAYLQYLYSEEGQRIAAQNFYRPRNQKVAAEFATQFPKLNLVTVDSDFGGWKTAQPKFFNDGGIFDQIYQAQ*</td><td>-</td><td>-</td><td>A0A071L1U4</td><td>⋯</td><td>99.698795</td><td>-1.000000</td><td>K02048: sbp; sulfate-binding protein                                      </td><td>Pseudomonas aeruginosa PAO1</td><td>PF12727.7; PF01547.25; PF13531.6; PF13343.6                                                                       </td><td>0.000540; 0.000000; 0.000000; 0.000000                                                            </td><td>PBP_like; SBP_bac_1; SBP_bac_11; SBP_bac_6                                      </td><td>TIGR00971                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         </td><td>0.000000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>TIGR00971                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         </td></tr>
	<tr><th scope=row>3619</th><td>caef037ad318582c885d129aee6008b1_3762</td><td>65.93</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MSRRISPVIPGFGLTLGYTLVYLSLLVLIPLGAMFLKTTQLSWEQFWAIISAPRVLAALKLSFGTALAAAVLNGLIGTLLAWVLVRYEFPGRKIIDAMIDLPFALPTAVAGIALTALYAPAGLVGQFASDLGFKIAYTPLGITLALTFVTLPFVVRTVQPVLADIPKEVEEAAACLGARPLQVFRHILVPALLPAWLTGFALAFARGVGEYGSVIFIAGNMPMKTEILPLLIMVKLDQYDYTGATAIGVLMLVVSFILLLLINLLQRRIETP*                                                            </td><td>-</td><td>-</td><td>A0A072ZRJ3</td><td>⋯</td><td>99.632353</td><td>-1.000000</td><td>K02046: cysT; sulfate transporter CysT                                    </td><td>Pseudomonas aeruginosa PAO1</td><td>PF00528.22                                                                                                        </td><td>0.000000                                                                                          </td><td>BPD_transp_1                                                                    </td><td>TIGR00969; TIGR00974; TIGR01253; TIGR01581; TIGR02138; TIGR02139; TIGR02140; TIGR02141; TIGR03226; TIGR03255; TIGR03262; TIGR03416                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000088                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>TIGR00969; TIGR00974; TIGR01253; TIGR01581; TIGR02138; TIGR02139; TIGR02140; TIGR02141; TIGR03226; TIGR03255; TIGR03262; TIGR03416                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td></tr>
	<tr><th scope=row>3618</th><td>caef037ad318582c885d129aee6008b1_3763</td><td>65.17</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MSSASISAATASASRRGNAGGRLALIVLAWLVFALFLLLPLYVVLSEALKQGFGTFFEAILEPDALAALKLTLIAVAISVPLNLLFGVAAAWCVSKFEFRGKSILVTLIDLPFSVSPVIAGLIYVLLFGAQGYFGEWLSDHDIQIVFAVPGIVLATLFVTVPFVARELIPLMQEQGTQEEEAARLLGANGWQMFWHVTLPNIKWGLIYGVVLCTARAMGEFGAVSVVSGHIRGVTNTLPLHVEILYNEYNHVAAFSVASLLLLMALVILLLKQWSESRMSRLKSNADEE*                                           </td><td>-</td><td>-</td><td>A0A069QBD9</td><td>⋯</td><td>99.653979</td><td>-1.000000</td><td>K02047: cysW; sulfate transporter CysW                                    </td><td>Pseudomonas aeruginosa PAO1</td><td>PF00528.22                                                                                                        </td><td>0.000000                                                                                          </td><td>BPD_transp_1                                                                    </td><td>TIGR00969; TIGR00974; TIGR01097; TIGR01253; TIGR01581; TIGR02138; TIGR02139; TIGR02140; TIGR02141; TIGR03226; TIGR03255; TIGR03262                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>0.000000; 0.000000; 0.000860; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>TIGR00969; TIGR00974; TIGR01097; TIGR01253; TIGR01581; TIGR02138; TIGR02139; TIGR02140; TIGR02141; TIGR03226; TIGR03255; TIGR03262                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td></tr>
	<tr><th scope=row>3625</th><td>caef037ad318582c885d129aee6008b1_3764</td><td>65.66</td><td>caef037ad318582c885d129aee6008b1</td><td>0</td><td>0</td><td>Reverse</td><td>MSIEIRNVSKNFNAFKALDNINLDIQSGELVALLGPSGCGKTTLLRIIAGLETPDAGNIVFHGEDVSQHDVRDRNVGFVFQHYALFRHMTVFDNVAFGLRMKPKGERPGESAIKAKVHELLNMVQLDWLADRYPEQLSGGQRQRIALARALAVEPKILLLDEPFGALDAKVRKELRRWLARLHEEINLTSVFVTHDQEEAMEVADRIVVMNKGVIEQIGSPGEVYENPASDFVYHFLGDSNRLQLGNDQHLLFRPHEVSLSRSAVAEHRAAEVRDIRPLGAITRVTLKVDGQDELIEAEVVKDHDSLAGLARGETLYFKPKAFQPVANL*   </td><td>-</td><td>-</td><td>A0A072ZFC4</td><td>⋯</td><td>99.696049</td><td>-1.000000</td><td>K02045: cysA; sulfate.thiosulfate ABC transporter ATP-binding protein CysA</td><td>Pseudomonas aeruginosa PAO1</td><td>PF00004.29; PF13191.6; PF13304.6; PF13476.6; PF00005.27; PF03215.15; PF03193.16; PF02463.19; PF08402.10; PF12857.7</td><td>0.000042; 0.000081; 0.000000; 0.000016; 0.000000; 0.000270; 0.000620; 0.000000; 0.000054; 0.000000</td><td>AAA; AAA_16; AAA_21; AAA_23; ABC_tran; Rad17; RsgA_GTPase; SMC_N; TOBE_2; TOBE_3</td><td>TIGR00602; TIGR00611; TIGR00618; TIGR00630; TIGR00954; TIGR00955; TIGR00956; TIGR00957; TIGR00958; TIGR00968; TIGR00972; TIGR01166; TIGR01184; TIGR01186; TIGR01187; TIGR01188; TIGR01189; TIGR01192; TIGR01193; TIGR01194; TIGR01257; TIGR01271; TIGR01277; TIGR01288; TIGR01842; TIGR01846; TIGR01978; TIGR02142; TIGR02203; TIGR02204; TIGR02211; TIGR02314; TIGR02315; TIGR02323; TIGR02324; TIGR02633; TIGR02673; TIGR02769; TIGR02770; TIGR02857; TIGR02868; TIGR02982; TIGR03005; TIGR03258; TIGR03265; TIGR03269; TIGR03375; TIGR03410; TIGR03411; TIGR03415; TIGR03522; TIGR03608; TIGR03719; TIGR03740; TIGR03771; TIGR03796; TIGR03797; TIGR03864; TIGR03873; TIGR04406</td><td>0.000930; 0.000380; 0.000760; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000; 0.000000</td><td>TIGR00602; TIGR00611; TIGR00618; TIGR00630; TIGR00954; TIGR00955; TIGR00956; TIGR00957; TIGR00958; TIGR00968; TIGR00972; TIGR01166; TIGR01184; TIGR01186; TIGR01187; TIGR01188; TIGR01189; TIGR01192; TIGR01193; TIGR01194; TIGR01257; TIGR01271; TIGR01277; TIGR01288; TIGR01842; TIGR01846; TIGR01978; TIGR02142; TIGR02203; TIGR02204; TIGR02211; TIGR02314; TIGR02315; TIGR02323; TIGR02324; TIGR02633; TIGR02673; TIGR02769; TIGR02770; TIGR02857; TIGR02868; TIGR02982; TIGR03005; TIGR03258; TIGR03265; TIGR03269; TIGR03375; TIGR03410; TIGR03411; TIGR03415; TIGR03522; TIGR03608; TIGR03719; TIGR03740; TIGR03771; TIGR03796; TIGR03797; TIGR03864; TIGR03873; TIGR04406</td></tr>
</tbody>
</table>

We found 2 genes that match to K02048: the sulfate binding protein. However, as we are only interested in the gene cluster, we will look at the one from the same contig as the other genes (caef037ad318582c885d129aee6008b1). 

Next, we will create a table (`*cys.txt`) to be used for `genoPlotR`.

```R
GOI2=GOI[2:5,] %>%  
    select("Description.2","Query.Gene") %>% 
    separate(col = Description.2, into = c("KO", "Annotation")) 

rownames(GOI2) <- c() 

head(GOI2)
```

    Warning message:
    “Expected 2 pieces. Additional pieces discarded in 4 rows [1, 2, 3, 4].”


<table>
<caption>A data.frame: 4 × 3</caption>
<thead>
	<tr><th scope=col>KO</th><th scope=col>Annotation</th><th scope=col>Query.Gene</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>K02048</td><td>sbp </td><td>caef037ad318582c885d129aee6008b1_3761</td></tr>
	<tr><td>K02046</td><td>cysT</td><td>caef037ad318582c885d129aee6008b1_3762</td></tr>
	<tr><td>K02047</td><td>cysW</td><td>caef037ad318582c885d129aee6008b1_3763</td></tr>
	<tr><td>K02045</td><td>cysA</td><td>caef037ad318582c885d129aee6008b1_3764</td></tr>
</tbody>
</table>

Export the file. 

```R
write_delim(GOI2, "bin2_cys.txt", delim="\t")
```