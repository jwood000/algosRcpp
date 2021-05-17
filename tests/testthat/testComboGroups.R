context("testing comboGroup")

test_that("comboGroups produces correct results", {
    expect_equal(nrow(comboGroups(letters[1:12], 4, "3Darray")), comboGroupsCount(12, 4))
    expect_equal(nrow(comboGroups(factor(letters[1:12]), 3, "3Darray")),
                 comboGroupsCount(12, 3))
    expect_equal(nrow(comboGroups(12, 1, "3Darray")), comboGroupsCount(12, 1))
    expect_equal(nrow(comboGroups(12, 12, "3Darray")), comboGroupsCount(12, 12))
    
    # comboGroupsCount(50, 5)
    # Big Integer ('bigz') :
    # [1] 402789797982510165934296910320
    expect_equal(dim(comboGroups(50, 5, "3Darray", lower = "402789797982510165934296910301")),
                 c(20, 10, 5))
    expect_equal(comboGroups(10, 5, lower = 201, upper = 300),
                 comboGroups(10, 5)[201:300, ])
    
    expect_equal(comboGroups(12, 4), 
                 comboGroupsSample(12, 4, sampleVec = 1:comboGroupsCount(12, 4)))
    
    expect_equal(comboGroups(12, 4), 
                 comboGroups(12, 4, nThreads = 2))
    expect_equal(comboGroups(12, 4), 
                 comboGroupsSample(12, 4, sampleVec = 1:comboGroupsCount(12, 4), nThreads = 2))
    
    expect_equal(dim(comboGroups(200, 5, "3Darray", upper = 100)),
                 c(100, 40, 5))
    
    expect_equal(comboGroups(14, 7, nThreads = 2), comboGroups(14, 7))
    expect_equal(sum(comboGroups(c(T, F), 2, "3Darray")), 1)
    
    expect_equal(rownames(comboGroupsSample(9, 3, "3Darray",
                                            sampleVec = c(67, 15, 248),
                                            namedSample = TRUE)),
                 as.character(c(67, 15, 248)))
    
    expect_equal(dim(comboGroups(as.raw(1:4), 2, "3Darray")),
                 dim(comboGroups(as.complex(c(1, -1, 1i, -1i)), 2, "3Darray")))
    
    expect_equal(rownames(comboGroupsSample(30, 5, n = 2,
                                            seed = 1, namedSample = TRUE)),
                 c("7162662695786451", "3525427663529072"))
    
    expect_equal(comboGroups(1:4 + 0.1, 2, "3Darray"), 
                 comboGroupsSample(1:4 + 0.1, 2, "3Darray", sampleVec = 1:3))
})

