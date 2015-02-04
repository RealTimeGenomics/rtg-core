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
package com.rtg.sam;

import java.util.Arrays;

import com.rtg.alignment.CgGotohEditDistance;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.util.Utils;

import net.sf.samtools.SAMRecord;

/**
 * Takes care of some of the legacy Cg handling for the SamValidator
 */
public final class SamValidatorCgHelper {

  private SamValidatorCgHelper() { }

  static byte[] expandCgCigarQualities(byte[] samQualities, byte[] qualityBuffer, String gq, int[] cggc, boolean first, boolean reverseCompliment, boolean phredifyQualityOverlap) {
    assert qualityBuffer.length == SamUtils.CG_RAW_READ_LENGTH;
    final int overlapQualityOffset = phredifyQualityOverlap ? FastaUtils.PHRED_LOWER_LIMIT_CHAR : 0;
    final int xqlength = gq == null ? 0 : gq.length();
    if (samQualities.length + xqlength != qualityBuffer.length) {
      return null;
    }
    if (gq == null || gq.length() == 0) {
      assert samQualities.length == qualityBuffer.length;
      System.arraycopy(samQualities, 0, qualityBuffer, 0, samQualities.length);
    } else {
      if ((first && !reverseCompliment) || (reverseCompliment && !first)) { //overlap on left hand side
        System.arraycopy(samQualities, 0, qualityBuffer, 0, cggc[0]);
        for (int i = 0; i < cggc[1]; i++) {
          qualityBuffer[cggc[0] + i] = (byte) (gq.charAt(i) - overlapQualityOffset);
        }
        System.arraycopy(samQualities, cggc[0], qualityBuffer, cggc[0] + cggc[1], cggc[2] + cggc[1]);
      } else {  //overlap on right hand side
        System.arraycopy(samQualities, 0, qualityBuffer, 0, cggc[0] + cggc[1]);
        for (int i = 0; i < cggc[1]; i++) {
          qualityBuffer[cggc[0] + cggc[1] + i] = (byte) (gq.charAt(i) - overlapQualityOffset);
        }
        System.arraycopy(samQualities, cggc[0] + cggc[1], qualityBuffer, cggc[0] + cggc[1] + cggc[1], cggc[2]);
      }
    }

    if (reverseCompliment) {
      Utils.reverseInPlace(qualityBuffer);
    }
    return qualityBuffer;
  }

  private static byte[] sQualityBuf = new byte[35];

  static boolean matchesCg(byte[] read, byte[] quality, SAMRecord record, boolean ignoreFragmentSize) {

    final byte[] recordBytes = record.getReadBases();
    final int[] cggc = parseCGGCAttribute(record);
    if (cggc != null) {
      //check if GC field values are valid (see bug 1144)
      if (recordBytes.length != (cggc[0] + cggc[1] + cggc[2])) {
        return false;
      }
      if (CgGotohEditDistance.CG_RAW_READ_LENGTH != (cggc[0] + 2 * cggc[1] + cggc[2])) {
        return false;
      }
      if (!ignoreFragmentSize) {
        final boolean forward = record.getFirstOfPairFlag() ^ record.getReadNegativeStrandFlag();
        if (forward) {
          if (cggc[1] + cggc[0] != CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE) {
            return false;
          }
        } else {
          if (cggc[1] + cggc[2] != CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE) {
            return false;
          }
        }
      }
      final String gq = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY);
      final int gqlen = gq != null ? gq.length() : 0;
      if (quality != null) {
        if (quality.length != record.getBaseQualities().length + gqlen) {
          return false;
        } else if (gq != null && gqlen != cggc[1]) {
          return false;
        }
      }
    } else {
      if (read.length != recordBytes.length /*+ CgGotohEditDistance.CG_INSERT_REGION_DEFAULT_SIZE*/) { //the read hasn't been aligned so will only contain the default number of cg gaps
        return false;
      } else if (quality != null && quality.length != record.getBaseQualities().length) {
        return false;
      }
    }
    final int cgOffset = 0;
    if (record.getReadNegativeStrandFlag()) {
      int cgOverlapSize = 0;
      for (int i = 0; i + cgOverlapSize < read.length; i++) {
//        if (read[i + cgOverlapSize] == DnaUtils.CG_SPACER_VALUE) {
//          cgOffset++;
//          continue;
//        }
        if (cggc != null && i - cgOffset == cggc[2]) {
          //we're at the start of the overlap region, determine if the ENTIRE overlap region is correct
          final String gs = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
          if (gs == null || gs.length() != 2 * cggc[1]) {
            return false;
          }
          if (!checkCGOverlapRegion(record, new int[] {cgOffset, cgOverlapSize}, i, read, quality, cggc)) {
            return false;
          }

          cgOverlapSize = gs.length() / 2; //skip that number of places along the read so we can continue as normal.
        } else {
          if (DnaUtils.getBase(DNA.complement(read[i + cgOverlapSize])) != Character.toUpperCase((char) recordBytes[recordBytes.length - i - 1 + cgOffset])) {
            return false;
//          } else if (quality != null && i < record.getBaseQualities().length && quality[i + cgOverlapSize - cgOffset] != record.getBaseQualities()[recordBytes.length - i - 1 + cgOffset]) {
//            System.err.println(cgOverlapSize + " :  " + cgOffset);
//            System.err.println(i + " :: " + record.getBaseQualities().length + " :: " + quality[i + cgOverlapSize - cgOffset] + " :: " + record.getBaseQualities()[recordBytes.length - i - 1 + cgOffset]);
//            System.err.println("bad4");
//            return false;
          }
        }
      }
    } else {
      int cgOverlapSize = 0;
      for (int i = 0; i + cgOverlapSize < read.length; i++) {
//        if (read[i + cgOverlapSize] == DnaUtils.CG_SPACER_VALUE) {
//          cgOffset++;
//          continue;
//        }
        if (cggc != null && i - cgOffset == cggc[0]) {
          //we're at the start of the overlap region, determine if the ENTIRE overlap region is correct
          final String gs = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
          if (gs == null || gs.length() != 2 * cggc[1]) {
            return false;
          }
          cgOverlapSize = gs.length() / 2; //cgOverlapSize will now skip that number of places along the read so we can continue as normal.
          if (!checkCGOverlapRegion(record, new int[] {cgOffset, cgOverlapSize}, i, read, quality, cggc)) {
            return false;
          }
        } else if (DnaUtils.getBase(read[i + cgOverlapSize]) != Character.toUpperCase((char) recordBytes[i - cgOffset])) {
          return false;
//        } else if (quality != null && i < record.getBaseQualities().length && quality[i + cgOverlapSize - cgOffset] != record.getBaseQualities()[i - cgOffset]) {
//          return false;
        }
      }
    }

    // check qualities...
    if (quality != null && record.getBaseQualities() != null) {
      final byte[] qualExpanded = expandCgCigarQualities(record.getBaseQualities(), sQualityBuf, record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY), cggc, record.getFirstOfPairFlag(), record.getReadNegativeStrandFlag(), true);
      if (qualExpanded == null) {
        return false;
      }
      if (!Arrays.equals(quality, qualExpanded)) {
        return false;
      }
    }

    return true;
  }

  private static boolean checkCGOverlapRegion(SAMRecord record, int[] cgData, int readPos, byte[] read, byte[] qual, int[] cggc) {
    final int cgOffset = cgData[0];
    final int cgOverlapSize = cgData[1];
    final String gs = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
    final String gq = record.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY);

    if (record.getFirstOfPairFlag()) {
      for (int j = 0; j < gs.length(); j++) {
        if (Character.toUpperCase(gs.charAt(record.getReadNegativeStrandFlag() ? gs.length() - 1 - j : j)) != DnaUtils.getBase(record.getReadNegativeStrandFlag() ? DNA.complement(read[readPos + j]) : read[readPos + j])) {
          return false;
        }
        final int recordCharPos = record.getReadNegativeStrandFlag() ? record.getReadBases().length - 1 - readPos + cgOffset - (j - gs.length() / 2) : readPos - cgOffset + (j - cgOverlapSize);
        if (j >= cggc[1]) {

          if (Character.toUpperCase((char) record.getReadBases()[recordCharPos]) != Character.toUpperCase(gs.charAt(record.getReadNegativeStrandFlag() ? gs.length() - 1 - j : j))) {
            return false;
          }
          if (qual != null) {
            if ((record.getBaseQualities()[recordCharPos]) != qual[readPos + j]) {
              return false;
            }
          }
        } else {
          if (gq != null && qual != null) {
            if ((gq.getBytes()[record.getReadNegativeStrandFlag() ? gq.length() - 1 - j : j]) - 33 != qual[readPos + j]) {
              return false;
            }
          }
        }
      }
    } else {
      final int recordCharPos = record.getReadNegativeStrandFlag() ? record.getReadBases().length - readPos - 1 + cgOffset : readPos - cgOffset;
      for (int j = 0; j < gs.length(); j++) {
        if (Character.toUpperCase(gs.charAt(record.getReadNegativeStrandFlag() ? gs.length() - 1 - j : j)) != DnaUtils.getBase(record.getReadNegativeStrandFlag() ? DNA.complement(read[readPos + j]) : read[readPos + j])) {
          return false;
        }
        if (j < cggc[1]) {
          if (DnaUtils.getBase(record.getReadNegativeStrandFlag() ? DNA.complement(read[readPos]) : read[readPos]) != Character.toUpperCase((char) record.getReadBases()[recordCharPos])) {
            return false;
          }
          if (qual != null) {
            if ((record.getBaseQualities()[recordCharPos - (record.getReadNegativeStrandFlag() ? j : -j)]) != qual[readPos + j - cgOffset]) {
              return false;
            }
          }
        } else {  //check the read quality vs the gq value...
          if (gq != null && qual != null) {
            if ((gq.getBytes()[record.getReadNegativeStrandFlag() ? gq.length() - 1 - (j - cggc[1]) : j - cggc[1]]) - 33 != qual[readPos + j - cgOffset]) {
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  static int[] parseCGGCAttribute(final SAMRecord samRecord) {
    final String cggc = samRecord.getStringAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
    if (cggc == null) {
      return null;
    }
    final int[] data = new int[3];
    int n = 0;

    for (int i = 0; i < cggc.length(); i++) {
      final char c = cggc.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        assert n > 0;
        switch (c) {
          case 'S':
            if (data[0] == 0) {
              data[0] = n;
            } else {
              data[2] = n;
            }
            break;
          case 'G':
            data[1] = n;
            break;
          default:
            assert false;
        }
        n = 0;
      }
    }
    assert data[0] > 0;
    assert data[1] > 0;
    assert data[2] > 0;
    return data;
  }
}
