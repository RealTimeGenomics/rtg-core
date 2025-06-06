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
package com.rtg.ngs;


import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.Arm;
import com.rtg.util.License;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Keeps track of information for determination of unmapped etymology
 */
public class ReadStatusTracker implements ReadStatusListener {

  /** Left arm is matched  */
  public static final int MATCHED_FIRST = 0x01;
  /** Right arm is matched */
  public static final int MATCHED_SECOND = 0x02;
  /** Left arm is blocked by read freq  */
  public static final int BLOCKED_FIRST = 0x04;
  /** Right arm is blocked by read freq  */
  public static final int BLOCKED_SECOND = 0x08;
  /** Both arms are matched and read is mated */
  public static final int MATED = 0x10;
  /** Read has been mated and <code>combo</code> score has been computed */
  public static final int MATED_ALIGN_SCORE = 0x20;
  /** Unmated left arm has alignment computed  */
  public static final int UNMATED_COMPUTE_ALIGNMENT_FIRST = 0x40;
  /** Unmated right arm has alignment computed  */
  public static final int UNMATED_COMPUTE_ALIGNMENT_SECOND = 0x80;
  /** Unmated left arm has passed <code>topn</code> threshold  */
  public static final int UNMATED_ALIGN_SCORE_FIRST = 0x100;
  /** Unmated right arm has passed <code>topn</code> threshold  */
  public static final int UNMATED_ALIGN_SCORE_SECOND = 0x200;
  /** First arm of read was written mated */
  public static final int MATED_FIRST = 0x400;
  /** Second arm of read was written mated */
  public static final int MATED_SECOND = 0x800;
  /** First arm of read was written unmated */
  public static final int UNMATED_FIRST = 0x1000;
  /** Second arm of read was written unmated */
  public static final int UNMATED_SECOND = 0x2000;
  /** First arm of read was written unmapped */
  public static final int UNMAPPED_FIRST = 0x4000;
  /** Second arm of read was written unmapped */
  public static final int UNMAPPED_SECOND = 0x8000;
  /** Single-end read was written unmapped */
  public static final int UNMAPPED = UNMAPPED_FIRST | UNMAPPED_SECOND;
  /** Uniquely mapped first arm */
  public static final int UNIQUELY_MAPPED_FIRST = 0x10000;
  /** Uniquely mapped second arm */
  public static final int UNIQUELY_MAPPED_SECOND = 0x20000;
  /** Zero-length (or shorter than threshold) first arm */
  public static final int SHORT_FIRST = 0x40000;
  /** Zero-length (or shorter than threshold) second arm */
  public static final int SHORT_SECOND = 0x80000;


  static final int MAPPING_STATUS_MASK = UNMATED_FIRST | UNMATED_SECOND | MATED_FIRST | MATED_SECOND | UNMAPPED_FIRST | UNMAPPED_SECOND;

  protected final MapStatistics mStatistics;

  protected final int[] mReadIdStatus;


  protected ReadStatusTracker(int numReads, MapStatistics stats) {
    mReadIdStatus = new int[numReads];
    mStatistics = stats;
  }

  /**
   * Add status to read
   * @param readId read id
   * @param status status to set
   */
  @Override
  public void addStatus(int readId, int status) {
    if (isSet(mReadIdStatus[readId], status)) {
      return;
    }
    mReadIdStatus[readId] |= status;
  }

  /**
   * Get the status of a particular read
   * @param readId read id
   * @param status status to check
   * @return true if given status is set for the given read id
   */
  public boolean getStatus(int readId, int status) {
    return isSet(mReadIdStatus[readId], status);
  }

  /**
   * Get a string representation of a status.
   * @param st status to describe
   * @return a string representation of the status
   */
  public static String statusToString(int st) {
    final StringBuilder str = new StringBuilder();
    if ((st & MATED_FIRST) != 0) {
      str.append("MATED_FIRST ");
    }
    if ((st & MATED_SECOND) != 0) {
      str.append("MATED_SECOND ");
    }
    if ((st & UNMATED_FIRST) != 0) {
      str.append("UNMATED_FIRST ");
    }
    if ((st & UNMATED_SECOND) != 0) {
      str.append("UNMATED_SECOND ");
    }
    if ((st & UNMAPPED_FIRST) != 0) {
      str.append("UNMAPPED_FIRST ");
    }
    if ((st & UNMAPPED_SECOND) != 0) {
      str.append("UNMAPPED_SECOND ");
    }
    if ((st & BLOCKED_FIRST) != 0) {
      str.append("BLOCKED_FIRST ");
    }
    if ((st & BLOCKED_SECOND) != 0) {
      str.append("BLOCKED_SECOND ");
    }
    if ((st & MATCHED_FIRST) != 0) {
      str.append("MATCHED_FIRST ");
    }
    if ((st & MATCHED_SECOND) != 0) {
      str.append("MATCHED_SECOND ");
    }
    if ((st & MATED_ALIGN_SCORE) != 0) {
      str.append("MATED_ALIGN_SCORE ");
    }
    if ((st & MATED) != 0) {
      str.append("MATED ");
    }
    if ((st & UNIQUELY_MAPPED_FIRST) != 0) {
      str.append("UNIQUELY_MAPPED_FIRST ");
    }
    if ((st & UNIQUELY_MAPPED_SECOND) != 0) {
      str.append("UNIQUELY_MAPPED_SECOND ");
    }
    if ((st & UNMATED_ALIGN_SCORE_FIRST) != 0) {
      str.append("UNMATED_ALIGN_SCORE_FIRST ");
    }
    if ((st & UNMATED_ALIGN_SCORE_SECOND) != 0) {
      str.append("UNMATED_ALIGN_SCORE_SECOND ");
    }
    if ((st & UNMATED_COMPUTE_ALIGNMENT_FIRST) != 0) {
      str.append("UNMATED_COMPUTE_ALIGNMENT_FIRST ");
    }
    if ((st & UNMATED_COMPUTE_ALIGNMENT_SECOND) != 0) {
      str.append("UNMATED_COMPUTE_ALIGNMENT_SECOND ");
    }
    if ((st & SHORT_FIRST) != 0) {
      str.append("SHORT_FIRST ");
    }
    if ((st & SHORT_SECOND) != 0) {
      str.append("SHORT_SECOND ");
    }
    return str.toString();
  }

  protected void calculateStatistics(boolean pairedEnd, boolean allhits) {
    if (mStatistics != null) {
      for (int r = 0; r < mReadIdStatus.length; ++r) {
        final int mappingStatus = mReadIdStatus[r] & MAPPING_STATUS_MASK;
        final int status = mReadIdStatus[r];
        boolean leftNoHit = false;
        if (isSet(status, SHORT_FIRST)) {
          mStatistics.increment(MapStatisticsField.IGNORED, Arm.LEFT);
        } else {
          mStatistics.increment(MapStatisticsField.TOTAL_READS, Arm.LEFT);
          if (isSet(status, MATED_FIRST)) {
            if (isSet(status, UNIQUELY_MAPPED_FIRST)) {
              mStatistics.increment(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
            } else {
              mStatistics.increment(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT);
            }
          } else if (isSet(status, UNMATED_FIRST)) {
            if (isSet(status, UNIQUELY_MAPPED_FIRST)) {
              mStatistics.increment(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT);
            } else {
              mStatistics.increment(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT);
            }
          } else {
            leftNoHit = true;
            if (isSet(status, BLOCKED_FIRST)) {
              mStatistics.increment(MapStatisticsField.UNMAPPED_BLOCKED, Arm.LEFT);
            } else if (isSet(status, MATED) && !(isSet(status, BLOCKED_FIRST) || isSet(status, BLOCKED_SECOND))) {
              if (!isSet(status, MATED_ALIGN_SCORE)) {
                mStatistics.increment(MapStatisticsField.UNMAPPED_MATED_POOR, Arm.LEFT);
              } else {
                mStatistics.increment(MapStatisticsField.UNMAPPED_MATED_TOO_MANY, Arm.LEFT);
              }
            } else {
              if (allhits) {
                if (isSet(status, MATCHED_FIRST) && isSet(status, UNMATED_COMPUTE_ALIGNMENT_FIRST)) {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT);
                } else {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
                }
              } else if (isSet(status, MATCHED_FIRST)) {
                if (!isSet(status, UNMATED_COMPUTE_ALIGNMENT_FIRST)) {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_TOPN, Arm.LEFT);
                } else if (!isSet(status, UNMATED_ALIGN_SCORE_FIRST)) {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT);
                } else {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, Arm.LEFT);
                }
              } else {
                mStatistics.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
              }
            }
          }
        }
        if (pairedEnd) {
          if (isSet(status, SHORT_SECOND)) {
            mStatistics.increment(MapStatisticsField.IGNORED, Arm.RIGHT);
          } else {
            mStatistics.increment(MapStatisticsField.TOTAL_READS, Arm.RIGHT);
            if (isSet(status, MATED_SECOND)) {
              if (isSet(status, UNIQUELY_MAPPED_SECOND)) {
                mStatistics.increment(MapStatisticsField.MATED_UNIQUE_READS, Arm.RIGHT);
              } else {
                mStatistics.increment(MapStatisticsField.MATED_AMBIG_READS, Arm.RIGHT);
              }
            } else if (isSet(status, UNMATED_SECOND)) {
              if (isSet(status, UNIQUELY_MAPPED_SECOND)) {
                mStatistics.increment(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.RIGHT);
              } else {
                mStatistics.increment(MapStatisticsField.UNMATED_AMBIG_READS, Arm.RIGHT);
              }
            } else {
              if (leftNoHit) {
                ((PairedEndMapStatistics) mStatistics).incrementBothUnmapped();
              }
              if (isSet(status, BLOCKED_SECOND)) {
                mStatistics.increment(MapStatisticsField.UNMAPPED_BLOCKED, Arm.RIGHT);
              } else if (isSet(status, MATED) && !(isSet(status, BLOCKED_FIRST) || isSet(status, BLOCKED_SECOND))) {
                if (!isSet(status, MATED_ALIGN_SCORE)) {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_MATED_POOR, Arm.RIGHT);
                } else {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_MATED_TOO_MANY, Arm.RIGHT);
                }
              } else {
                if (allhits) {
                  if (isSet(status, MATCHED_SECOND) && isSet(status, UNMATED_COMPUTE_ALIGNMENT_SECOND)) {
                    mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.RIGHT);
                  } else {
                    mStatistics.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.RIGHT);
                  }
                } else if (isSet(status, MATCHED_SECOND)) {
                  if (!isSet(status, UNMATED_COMPUTE_ALIGNMENT_SECOND)) {
                    mStatistics.increment(MapStatisticsField.UNMAPPED_TOPN, Arm.RIGHT);
                  } else if (!isSet(status, UNMATED_ALIGN_SCORE_SECOND)) {
                    mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.RIGHT);
                  } else {
                    mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, Arm.RIGHT);
                  }
                } else {
                  mStatistics.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.RIGHT);
                }
              }
            }
          }
        }

        if (License.isDeveloper()) {
          switch (mappingStatus) {
          case MATED_FIRST | MATED_SECOND:
          case UNMAPPED_SECOND | UNMATED_FIRST: //
          case UNMAPPED_FIRST | UNMATED_SECOND: //
          case UNMATED_FIRST | UNMATED_SECOND: //
          case UNMATED_FIRST: // should only happen when outputting unmated while not outputting unmapped
          case UNMATED_SECOND: //
          case UNMAPPED_FIRST | UNMAPPED_SECOND: //writing unmapped
          case 0x00: //is unmapped but we're not recording unmapped
            break;
          case MATED_FIRST: //missing right arm
          case UNMAPPED_FIRST: //
            mStatistics.increment(MapStatisticsField.MISSING, Arm.LEFT);
            Diagnostic.developerLog("readId: " + r + " had no right side status. Status was: " + statusToString(mappingStatus));
            break;
            // missing left arm
          case MATED_SECOND: //
          case UNMAPPED_SECOND: //
            mStatistics.increment(MapStatisticsField.MISSING, Arm.RIGHT);
            Diagnostic.developerLog("readId: " + r + " had no left side status. Status was: " + statusToString(mappingStatus));
            break;
          default:
            //error
            Diagnostic.developerLog("readId: " + r + " has faulty status: " + statusToString(mappingStatus) + "\n");
            break;
          }
        }
      }
    }
    Diagnostic.userLog("Setting stats: " + mStatistics);
  }

  boolean isSet(int val, int attr) {
    return (val & attr) == attr;
  }

  /**
   * Returns the appropriate <code>XC</code> attribute for read
   * @param readId the ID of the read
   * @param first whether is left read or not
   * @return <code>XC code</code> or <code>'\0'</code> if none
   */
  public char getXCAttribute(int readId, boolean first) {
    final int val = mReadIdStatus[readId];
    if ((first && isSet(val, BLOCKED_FIRST)) || (!first && isSet(val, BLOCKED_SECOND))) {
      return 'B';
    } else if (isSet(val, MATED) && !(isSet(val, BLOCKED_FIRST) || isSet(val, BLOCKED_SECOND))) {
      if (!isSet(val, MATED_ALIGN_SCORE)) {
        return 'd';
      }
      return 'e';
    } else {
      final int matchedAttr;
      final int computeAttr;
      final int alignScoreAttr;
      if (first) {
        matchedAttr = MATCHED_FIRST;
        computeAttr = UNMATED_COMPUTE_ALIGNMENT_FIRST;
        alignScoreAttr = UNMATED_ALIGN_SCORE_FIRST;
      } else {
        matchedAttr = MATCHED_SECOND;
        computeAttr = UNMATED_COMPUTE_ALIGNMENT_SECOND;
        alignScoreAttr = UNMATED_ALIGN_SCORE_SECOND;
      }
      if (isSet(val, matchedAttr)) {
        if (!isSet(val, computeAttr)) {
          return 'C';
        }
        if (!isSet(val, alignScoreAttr)) {
          return 'D';
        }
        return 'E';
      }
    }
    return 'A';
    //return '\0';
  }

  /**
   * @return the number of reads
   */
  public int getNumReads() {
    return mReadIdStatus.length;
  }

  /**
   * Ensure status for unmapped reads is correctly set
   * @param paired true if this is paired-end data
   */
  void preProcessUnMappedStatistics(boolean paired) {
    Diagnostic.progress("UnmappedPreprocess: Starting 1 Jobs");
    for (int i = 0; i < mReadIdStatus.length; ++i) {
      if (paired) {
        final boolean leftUnmapped = (mReadIdStatus[i] & (MATED_FIRST | UNMATED_FIRST)) == 0;
        final boolean rightUnmapped = (mReadIdStatus[i] & (MATED_SECOND | UNMATED_SECOND)) == 0;
        if (leftUnmapped) {
          addStatus(i, UNMAPPED_FIRST);
        }
        if (rightUnmapped) {
          addStatus(i, UNMAPPED_SECOND);
        }
      } else {
        if ((mReadIdStatus[i] & (MATED_FIRST | UNMATED_FIRST)) == 0) {
          addStatus(i, UNMAPPED);
        }
      }
    }
    Diagnostic.progress("UnmappedPreprocess: 1/1 Jobs Finished");
  }

  UnmappedStatus getUnmappedStatus(int readId, boolean paired) {
    if (paired) {
      final boolean leftUnmapped = (mReadIdStatus[readId] & UNMAPPED_FIRST) != 0;
      final boolean rightUnmapped = (mReadIdStatus[readId] & UNMAPPED_SECOND) != 0;
      if (leftUnmapped && rightUnmapped) {
        return UnmappedStatus.BOTH_UNMAPPED;
      }
      if (leftUnmapped) {
        return UnmappedStatus.LEFT_UNMAPPED;
      }
      if (rightUnmapped) {
        return UnmappedStatus.RIGHT_UNMAPPED;
      }
    } else {
      if ((mReadIdStatus[readId] & UNMAPPED) != 0) {
        return UnmappedStatus.SINGLE_END_UNMAPPED;
      }
    }
    return UnmappedStatus.MAPPED;
  }

  enum UnmappedStatus {
    LEFT_UNMAPPED,
    RIGHT_UNMAPPED,
    BOTH_UNMAPPED,
    SINGLE_END_UNMAPPED,
    MAPPED
  }
}
