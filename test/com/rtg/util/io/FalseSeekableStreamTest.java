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

package com.rtg.util.io;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * test class
 */
public class FalseSeekableStreamTest extends TestCase {

  private InputStream getStream() {
    final byte[] b = "012345678901234567890123456789".getBytes();
    return new ByteArrayInputStream(b);
  }

  public void testGetPosition() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.getPosition();
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testLength() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.length();
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testSeek() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.seek(5);
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testRead0args() throws Exception {
    final InputStream basic = getStream();
    try {
      final FalseSeekableStream fss = new FalseSeekableStream(getStream());
      try {
        assertEquals(basic.read(), fss.read());

      } finally {
        fss.close();
      }
    } finally {
      basic.close();
    }
  }

  public void testAvailable() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.available();
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testClose() throws Exception {
    final InputStream basic = getStream();
    try {
      final FalseSeekableStream fss = new FalseSeekableStream(getStream());
      fss.close();
    } finally {
      basic.close();
    }
  }

  public void testMark() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.mark(500);
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testMarkSupported() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.markSupported();
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testReadbyteArr() throws Exception {
    final InputStream basic = getStream();
    try {
      final FalseSeekableStream fss = new FalseSeekableStream(getStream());
      try {
        byte[] first = new byte[10];
        byte[] second = new byte[10];
        int len1 = basic.read(first);
        int len2 = fss.read(second);
        assertEquals(len1, len2);
        assertTrue(Arrays.equals(first, second));
      } finally {
        fss.close();
      }
    } finally {
      basic.close();
    }
  }

  public void testRead3args() throws Exception {
    final InputStream basic = getStream();
    try {
      final FalseSeekableStream fss = new FalseSeekableStream(getStream());
      try {
        byte[] first = new byte[10];
        byte[] second = new byte[10];
        int len1 = basic.read(first, 2, 8);
        int len2 = fss.read(second, 2, 8);
        assertEquals(len1, len2);
        assertTrue(Arrays.equals(first, second));
      } finally {
        fss.close();
      }
    } finally {
      basic.close();
    }
  }

  public void testReset() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      fss.reset();
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

  public void testSkip() throws IOException {
    final FalseSeekableStream fss = new FalseSeekableStream(getStream());
    try {
      long len = fss.skip(3);
      assertEquals(3, (int) len);
      fail();
    } catch (UnsupportedOperationException e) {
    } finally {
      fss.close();
    }
  }

}
