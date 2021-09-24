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

package com.rtg.variant;

import java.io.ByteArrayInputStream;
import java.io.File;

import com.rtg.tabix.TabixIndexer;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;

/**
 */
public final class VariantTestUtils {
  private VariantTestUtils() { }

  /**
   * Create SAM file and asociated tabix index from a resource given a directory where to put results.
   * @param resource path to a resource file.
   * @param dir sam file.
   * @return the SAM file.
   * @throws Exception whenever.
   */
  public static File bgzipAndIndexResource(final String resource, final File dir) throws Exception {
    final String samString = FileHelper.resourceToString(resource);
    final File sam = new File(dir, "sam.sam.gz");
    VariantTestUtils.bgzipAndIndex(samString, sam);
    return sam;
  }

  /**
   * Create SAM file and asociated tabix index.
   * @param sam text of a SAM file.
   * @param out sam file.
   * @throws Exception whenever.
   */
  public static void bgzipAndIndex(final String sam, final File out) throws Exception {
    BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(sam.getBytes()), out);
    new TabixIndexer(out, new File(out.getParent(), out.getName() + ".tbi")).saveSamIndex();
  }

}
