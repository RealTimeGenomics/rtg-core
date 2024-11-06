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

package com.rtg.index;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.AbstractTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Test
 */
public class HashBlacklistTest extends AbstractTest {

  public void testBlacklistExists() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File blacklistDir = new File(dir, HashBlacklist.BLACKLIST_SUBDIR);
      assertTrue(blacklistDir.mkdir());
      final File w15File = new File(blacklistDir, "w15");
      assertFalse(HashBlacklist.blacklistExists(dir, 15));
      assertTrue(w15File.createNewFile());
      assertTrue(HashBlacklist.blacklistExists(dir, 15));
    }
  }

  private static final String EXAMPLE_BLACKLIST = ""
    + "AAAAAAAA\t80" + StringUtils.LS
    + "AGAGAGAG\t100" + StringUtils.LS
    + "TTTTTTTT\t500" + StringUtils.LS
    + "CTCTCTCT\t1000" + StringUtils.LS;

  private static final Long[] EXPECTED_HASHES = {
    0b0000000000000000L,
    0b0010001000100010L,
    0b1111111111111111L,
    0b0111011101110111L
  };

  public void testLoad() throws IOException {
    List<Long> res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 0);
    assertEquals(Arrays.asList(EXPECTED_HASHES), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 80);
    assertEquals(Arrays.asList(EXPECTED_HASHES), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 81);
    assertEquals(Arrays.asList(EXPECTED_HASHES[1], EXPECTED_HASHES[2], EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 101);
    assertEquals(Arrays.asList(EXPECTED_HASHES[2], EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 501);
    assertEquals(Collections.singletonList(EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 1001);
    assertEquals(Collections.emptyList(), res);
  }
}
