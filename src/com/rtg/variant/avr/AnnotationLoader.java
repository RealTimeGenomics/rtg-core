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

package com.rtg.variant.avr;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Loads and saves {@link Annotation}s.
 */
public final class AnnotationLoader {
  protected enum AnnotationType {
    QUAL {
      @Override
      Annotation loadAnnotation(DataInputStream dis) {
        return new QualAnnotation();
      }
    },
    FORMAT {
      @Override
      Annotation loadAnnotation(DataInputStream dis) throws IOException {
        return new FormatAnnotation(dis.readUTF(), AnnotationDataType.values()[dis.readInt()]);
      }
    },
    INFO {
      @Override
      Annotation loadAnnotation(DataInputStream dis) throws IOException {
        return new InfoAnnotation(dis.readUTF(), AnnotationDataType.values()[dis.readInt()]);
      }
    },
    DERIVED {
      @Override
      Annotation loadAnnotation(DataInputStream dis) throws IOException {
        return new DerivedAnnotation(dis.readUTF());
      }
    };

    /**
     * Load an annotation from the given data input stream.  Stream is expected to contain valid values.
     * @param dis data input stream
     * @return an annotation base on stream content
     * @throws IOException if an error occurs reading
     */
    abstract Annotation loadAnnotation(DataInputStream dis) throws IOException;
  }

  /**
   * Prevent instantiation.
   */
  private AnnotationLoader() {
  }

  /**
   * Load an annotation from the given data input stream.
   * @param dis data input stream
   * @return an annotation base on stream content
   * @throws IOException if an error occurs reading
   */
  public static Annotation load(DataInputStream dis) throws IOException {
    return AnnotationType.values()[dis.readInt()].loadAnnotation(dis);
  }
}
