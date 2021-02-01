# Testing the STK Assignment
devtools::load_all()

library(JustinKinomeModelling)

data_file <- "kinome_data/STK_test_data.txt"

stk_data <- PamchipData_STK(data_file)

ampk_reporter_peptides <- c(
  "GRIK2_708_720",
  "KAP2_92_104",
  "NCF1_321_333",
  "CFTR_730_742",
  "MYPC3_268_280",
  "CREB1_126_138",
  "LMNB1_16_28",
  "CFTR_761_773",
  "PTK6_436_448",
  "VTNC_390_402",
  "KPB1_1011_1023",
  "RS6_228_240",
  "DESP_2842_2854",
  "NCF1_296_308",
  "SCN7A_898_910"
)

filtered <- subset_data(stk_data, ampk_reporter_peptides, "CTL_UT")

mod <- make_model(filtered, threshold = 0.75)

assigned <- assign_kinases(mod$avg_net, "STK")
