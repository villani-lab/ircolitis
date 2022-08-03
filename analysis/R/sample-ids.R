# Dictionary to convert internal to public identifiers
x <- read.table(header = TRUE, text = "sic_id  pub_id
SIC_100 C1   
SIC_121 C2
SIC_126 C3
SIC_134 C4
SIC_140 C5
SIC_141 C6
SIC_32  C7
SIC_36  C8
SIC_40  C9
SIC_43  C10
SIC_71  C11
SIC_89  C12
SIC_97  C13
SIC_109 IHC1
SIC_172 IHC2
SIC_19  IHC3
SIC_31  IHC4
SIC_53  IHC5
SIC_94  IHC6
SIC_132 IHC7
SIC_33  IHC8
SIC_13  HC1
SIC_186 HC2
SIC_187 HC3
SIC_188 HC4
SIC_196 HC5
MC_1    HC6
MC_2    HC7
MC_9    HC8")
pub_ids <- x$pub_id
names(pub_ids) <- x$sic_id

