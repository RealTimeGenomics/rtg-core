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
