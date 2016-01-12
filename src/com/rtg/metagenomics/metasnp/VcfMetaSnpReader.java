/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.metagenomics.metasnp;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.util.io.FileUtils;
import com.rtg.vcf.VcfReader;

/**
 * Adapter to convert from VCF input into records for metasnp.
 */
class VcfMetaSnpReader extends VcfReader implements MetaSnpReader {

  VcfMetaSnpReader(final File in) throws IOException {
    this(new BufferedReader(FileUtils.createReader(in, false)));
  }

  VcfMetaSnpReader(final BufferedReader in) throws IOException {
    super(in);
  }

  @Override
  public List<String> samples() {
    return getHeader().getSampleNames();
  }

  @Override
  public MetaSnpLine nextLine() throws IOException {
    while (hasNext()) {
      final MetaSnpLine metaSnpLine = MetaSnpLine.create(next());
      if (metaSnpLine != null) {
        return metaSnpLine;
      }
    }
    return null;
  }
}
