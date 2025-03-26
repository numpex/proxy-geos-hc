;; Ce qui suit est un "manifeste" équivalent à la ligne de commande que vous avez donnée.
;; Vous pouvez le stocker dans un fichier que vous pourrez ensuite passer à n'importe quelle
;; commande 'guix' qui accepte une option '--manifest' (ou '-m').
(use-modules (guix transformations))

(define kokkostransform
  (options->transformation
    '((with-commit . "kokkos=4.0.01"))))

(concatenate-manifests (list (packages->manifest
  (list (kokkostransform
          (specification->package "kokkos-openmp")))) (specifications->manifest
  (list "clang-toolchain@14.0.6"
        "python-sphinx"
        "graphviz"
        "doxygen"
        "uncrustify"
        "git"
        "grep"
        "sed"
        "valgrind"
        "python-yapf"
        "cmakelang"
        "cppcheck"
        "make"
        "astyle"
        "lesspipe"
        "cuda-toolkit"
        "vim"
        "binutils"
        "coreutils"
        "plocate"
        "cmake@3.25"
        "gcc-toolchain@11"))))
