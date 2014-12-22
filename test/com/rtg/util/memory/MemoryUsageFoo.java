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
package com.rtg.util.memory;

/**
 * A Class to be used by the MemoryUsageTest. This was originally a
 * static inner class of the MemoryUsageTest but the Jode 1.1.1
 * obfuscator failed the test. The failure I suspect was due to not
 * handling inner classes properly. But moving the test class to a
 * proper class it is hoped the Jode obfuscator will work .
 *
 */

class MemoryUsageFoo {
  protected int mD1, mD2;
}
