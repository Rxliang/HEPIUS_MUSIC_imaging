function p = pdata2roifn(pd, lambda_mm)

p.xMin = pd.Origin(1);
p.zMin = pd.Origin(3);

p.xMin_mm = p.xMin*lambda_mm;
p.zMin_mm = p.zMin*lambda_mm;

p.xExt = pd.Size(2)*pd.PDelta(1);
p.zExt = pd.Size(1)*pd.PDelta(3);

p.xExt_mm = p.xExt*lambda_mm;
p.zExt_mm = p.zExt*lambda_mm;

p.xMax = p.xMin + p.xExt;
p.zMax = p.zMin + p.zExt;

p.xMax_mm = p.xMin_mm + p.xExt_mm;
p.zMax_mm = p.zMin_mm + p.zExt_mm;

