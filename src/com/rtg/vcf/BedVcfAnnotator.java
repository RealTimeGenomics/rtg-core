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

package com.rtg.vcf;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedRangeLoader;
import com.rtg.bed.BedRecord;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Annotates a VCF records that match BED regions with the text from the matching BED region.
 * Note: Only uses the values from the name BED column (column 4).
 */
@TestClass("com.rtg.vcf.VcfAnnotatorCliTest")
public class BedVcfAnnotator implements VcfAnnotator {

  private static final class AnnotatorBedLoader extends BedRangeLoader {
    private AnnotatorBedLoader() {
      super(1);
    }
    @Override
    public String getMeta(BedRecord rec) {
      return new String(rec.getAnnotations()[0].toCharArray());
    }
  }


  /** Per chromosome annotations. */
  private final ReferenceRanges mAnnotations;

  private final boolean mAddIds;

  /**
   * Constructor
   * @param addIds if true, BED annotations will be added to VCF id column.
   * @param bedFiles BED files containing annotations
   * @throws IOException if the BED file could not be loaded.
   */
  public BedVcfAnnotator(boolean addIds, Collection<File> bedFiles) throws IOException {
    mAddIds = addIds;
    final BedVcfAnnotator.AnnotatorBedLoader bedLoader = new AnnotatorBedLoader();
    bedLoader.loadRanges(bedFiles);
    mAnnotations = bedLoader.getReferenceRanges();
  }

  private List<String> annotate(String chr, int loc) {
    List<String> anno = null;
    final RangeList<String> chrRanges = mAnnotations.get(chr);
    if (chrRanges != null) {
      anno = chrRanges.find(loc);
    }
    return anno;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    if (!mAddIds) {
      header.ensureContains(new InfoField("ANN", MetaType.STRING, new VcfNumber("."), "Annotation"));
    }
  }

  @Override
  public void annotate(VcfRecord rec) {
    final List<String> annotation = annotate(rec.getSequenceName(), rec.getStart());
    if (annotation != null) {
      if (annotation.size() > 0) {
        if (mAddIds) {
          rec.setId(annotation.toArray(new String[annotation.size()]));
        } else {
          rec.addInfo("ANN", annotation.toArray(new String[annotation.size()]));
        }
      }
    }
  }
}
