data("convertMsp2Spectra", package = "MetCirc")

spl <- convertMsp2Spectra(msp = msp2spectra)
## START unit test convertMSP2MSP
test_that("convertMsp2Spectra", {
    expect_equal(length(spl), 22)
    expect_is(spl, "Spectra")
    expect_equal(length(spl@backend), 22)
    expect_true(is.character(spl$name))
    expect_equal(spl@metadata, list())
    expect_is(spl$intensity, "SimpleNumericList")
    expect_true(is.numeric(spl$intensity[[1]]))
    expect_equal(spl$intensity[[1]],
        c(19794, 990, 1477, 267, 292, 332, 16681, 772, 24215, 954, 361, 225, 
            835, 6216, 307, 834, 1032, 328, 313, 1380, 4097, 418, 4553, 273,
            1176, 232,289, 735, 1844, 622, 256, 5600, 429, 384, 683, 220, 1269,
            10060, 851, 6905, 522, 349, 1152, 312, 238, 65655, 5279, 234, 271,
            422))
    expect_equal(spl$intensity[[2]],
        c(2744, 2412, 3763, 2210, 2562, 2958, 1879391, 1879391, 9060, 14055,
            9813, 9347, 393317, 2171, 2103, 1985, 1983))
    expect_true(is.numeric(spl$intensity[[2]]))
    expect_is(spl$mz, "SimpleNumericList")
    expect_true(is.numeric(spl$mz[[1]]))
    expect_equal(spl$mz[[1]], 
        c(117.0362, 118.0391, 121.0309, 135.0465, 141.0356, 148.0162, 149.0250,
            150.0279, 151.0046, 152.0079, 155.0506, 156.0567, 157.0657,
            159.0457, 160.0494, 161.0248, 169.0515, 171.0447, 173.0611,
            180.0586, 181.0665, 182.0490, 183.0462, 184.0483, 185.0261,
            186.5647, 195.0458, 196.0536, 197.0607, 199.0402, 200.7310,
            201.0557, 202.0588, 209.0248, 213.0564, 223.0403, 224.0486,
            225.0560, 226.0535, 227.0349, 228.0390, 233.0910, 241.0496,
            251.0343, 268.0359, 269.0447, 270.0483, 274.5693, 275.0310,
            766.5048),
        tolerance = 1e-04)
    expect_true(is.numeric(spl$mz[[2]]))
    expect_equal(spl$mz[[2]], 
        c(150.955, 171.014, 210.889, 214.887, 221.139, 256.974, 267.985,
            267.989, 279.998, 280.485, 281.079, 281.262, 283.022, 358.111,
            381.231, 398.138, 467.196),
        tolerance = 1e-03)
    expect_true(is.numeric(spl$precursorMz))
    expect_equal(spl$precursorMz, 
        c(269.0455, 283.0612, 301.0718, 271.0612, 269.0819, 301.0354, 591.1719,
            315.0510, 431.0984, 433.0776, 269.0455, 253.0506, 269.0455,
            609.1825, 301.0718, 463.0882, 285.0405, 285.0405, 609.1461,
            445.1140, 593.1301, 593.1512), 
        tolerance = 1e-04)
})
## END unit test convertMSP2MSP
