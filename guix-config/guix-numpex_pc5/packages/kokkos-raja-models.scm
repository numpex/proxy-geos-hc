;;; This module extends GNU Guix and is licensed under the same terms, those
;;; of the GNU GPL version 3 or (at your option) any later version.
;;;
;;; Copyright © 2024 Inria

(define-module (guix-numpex_pc5 packages kokkos-raja-models)
  #:use-module ((guix licenses)
                #:prefix license:)
  #:use-module (guix gexp)
  #:use-module (guix utils)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system cmake)
  #:use-module (guix build-system gnu)
  #:use-module (gnu packages cpp)
  #:use-module (guix-hpc-non-free packages cpp)
  #:use-module (llnl tainted geos)
  #:use-module (guix-science-nonfree packages cuda))

;; This creates a kokkos version with specific tweaks for gyselalibxx / with OpenMP enabled
(define (make-kokkos-cuda-openmp name kokkos-cuda-arch)
  (package/inherit kokkos-cuda-arch
    (name name)
    (arguments (substitute-keyword-arguments (package-arguments kokkos-cuda-arch)
                 ((#:configure-flags flags)
                  #~(append (list "-DKokkos_ENABLE_OPENMP=ON")
                            #$flags))
                 ;; Cannot run tests due to lack of specific hardware
                 ((#:tests? _ #t)
                  #f)
                 ;; RUNPATH validation fails since libcuda.so.1 is not present at build
                 ;; time.
                 ((#:validate-runpath? #f #f)
                  #f)
                 ((#:phases phases
                   '%standard-phases)
                  #~(modify-phases #$phases
                      ;; File is not present in CUDA build
                      ))
		 ))
    ))

(define-public kokkos-cuda-k40-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-k40-openmp" kokkos-cuda-k40))

(define-public kokkos-cuda-a40-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-a40-openmp" kokkos-cuda-a40))

(define-public kokkos-cuda-a100-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-a100-openmp" kokkos-cuda-a100))

(define-public kokkos-cuda-v100-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-v100-openmp" kokkos-cuda-v100))

(define-public kokkos-cuda-p100-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-p100-openmp" kokkos-cuda-p100))

(define-public kokkos-cuda-ada-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-ada-openmp" kokkos-cuda-ada))

(define-public kokkos-cuda-t4-openmp
  (make-kokkos-cuda-openmp "kokkos-cuda-t4-openmp" kokkos-cuda-t4))


;; This creates a raja-cuda where openmp and cuda are enabled throughout inheritence from raja-cuda and specification of a different compute capability
(define (make-raja-cuda-spec-compute name cuda-arch-compute)
  (package/inherit raja-cuda
    (name name)
    (arguments (substitute-keyword-arguments (package-arguments raja-cuda)
                 ((#:configure-flags flags)
                  #~(append (list (string-append "-DCMAKE_CUDA_ARCHITECTURES="#$cuda-arch-compute))
			          (delete "-DCMAKE_CUDA_ARCHITECTURES=70" #$flags)
                            ))
		 ))
    ))
(define-public raja-cuda-ada
  (make-raja-cuda-spec-compute "raja-cuda-ada" "89"))
(define-public raja-cuda-v100
  (make-raja-cuda-spec-compute "raja-cuda-v100" "70"))
(define-public raja-cuda-t4
  (make-raja-cuda-spec-compute "raja-cuda-t4" "75"))
(define-public raja-cuda-p100
  (make-raja-cuda-spec-compute "raja-cuda-p100" "60"))
(define-public raja-cuda-k40
  (make-raja-cuda-spec-compute "raja-cuda-k40" "35"))
(define-public raja-cuda-a40
  (make-raja-cuda-spec-compute "raja-cuda-a40" "86"))
(define-public raja-cuda-a100
  (make-raja-cuda-spec-compute "raja-cuda-a100" "80"))
