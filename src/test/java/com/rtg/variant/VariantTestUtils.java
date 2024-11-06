/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
