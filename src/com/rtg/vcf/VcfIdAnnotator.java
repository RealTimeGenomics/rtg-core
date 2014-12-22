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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.RangeList.RangeData;
import com.rtg.vcf.header.VcfHeader;

/**
 * Adds VCF ID column based on VCF IDs from another VCF file.
 */
@TestClass("com.rtg.vcf.VcfAnnotatorCliTest")
public class VcfIdAnnotator implements VcfAnnotator {

  private final Map<String, RangeList<String>> mAnnotations;

  /**
   * Constructor
   * @param vcfFiles VCF files containing variant IDs to be added to VCF id column.
   * @throws IOException if the BED file could not be loaded.
   */
  public VcfIdAnnotator(Collection<File> vcfFiles) throws IOException {
    final Map<String, List<RangeData<String>>> annotations = new HashMap<>();

    for (final File vcfFile : vcfFiles) {
      try (final VcfReader reader = VcfReader.openVcfReader(vcfFile)) {
        while (reader.hasNext()) {
          final VcfRecord record = reader.next();
          final String id = record.getId();
          if (!VcfRecord.MISSING.equals(id)) {
            final String[] ids = StringUtils.split(id, ';');
            final List<RangeList.RangeData<String>> annos;
            if (annotations.containsKey(record.getSequenceName())) {
              annos = annotations.get(record.getSequenceName());
            } else {
              annos = new ArrayList<>();
              annotations.put(record.getSequenceName(), annos);
            }
            final int start = record.getStart();
            final int end = start + record.getRefCall().length();
            final RangeList.RangeData<String> range = new RangeData<>(start, end, ids[0]);
            annos.add(range);
            for (int i = 1; i < ids.length; i++) {
              range.addMeta(ids[i]);
            }
          }
        }
      }
    }
    final Map<String, RangeList<String>> annotationSearch = new HashMap<>();
    for (final Map.Entry<String, List<RangeList.RangeData<String>>> me : annotations.entrySet()) {
      annotationSearch.put(me.getKey(), new RangeList<>(me.getValue()));
    }
    mAnnotations = annotationSearch;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    //No header changes necessary, modifying VCF id column
  }

  private List<String> annotate(String chr, int loc) {
    List<String> anno = null;
    if (mAnnotations.containsKey(chr)) {
      anno = mAnnotations.get(chr).find(loc);
    }
    return anno;
  }

  @Override
  public void annotate(VcfRecord rec) {
    final List<String> annotation = annotate(rec.getSequenceName(), rec.getStart());
    if (annotation != null) {
      if (annotation.size() > 0) {
        rec.setId(annotation.toArray(new String[annotation.size()]));
      }
    }
  }

}
