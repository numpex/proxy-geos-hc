;; Ce qui suit est un "manifeste" équivalent à la ligne de commande que vous avez donnée.
;; Vous pouvez le stocker dans un fichier que vous pourrez ensuite passer à n'importe quelle
;; commande 'guix' qui accepte une option '--manifest' (ou '-m').

(use-modules (guix transformations))

(define transform1
  (options->transformation
    '((with-configure-flag
        .
        "kokkos-cuda-ada=-DKokkos_ENABLE_OPENMP=ON"))))

(packages->manifest
  (list (transform1
          (specification->package "kokkos-cuda-ada"))
        (specification->package "clang-toolchain@14")
        (specification->package "python-sphinx")
        (specification->package "graphviz")
        (specification->package "doxygen")
        (specification->package "uncrustify")
        (specification->package "git")
        (specification->package "grep")
        (specification->package "sed")
        (specification->package "valgrind")
        (specification->package "python-yapf")
        (specification->package "cmakelang")
        (specification->package "cppcheck")
        (specification->package "make")
        (specification->package "astyle")
        (specification->package "lesspipe")
        (specification->package "cuda-toolkit")
        (specification->package "vim")
        (specification->package "binutils")
        (specification->package "coreutils")
        (specification->package "plocate")
        (specification->package "cmake@3.25")
        (specification->package "gcc-toolchain@11")))
