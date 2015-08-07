
predict(mnom.length.bulk.order3.nl, newdata = data.frame(
                                   prop.Carapace.Length = cl, 
                                   prop.Carapace.Length2 = cl^2, 
                                   prop.Carapace.Length3 = cl^3,
                                   mesh.70mm_Total.Weight = 110,
                                   mesh.80mm_Total.Weight = 110,
                                   mesh.90mm_Total.Weight = 110,
                                   mesh.100mm_Total.Weight = 110
                                   ), type = "prob")
