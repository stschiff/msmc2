import model.gsl_matrix_vector;
import std.stdio;

class SubVec {
    gsl_vector_view _view;
    this(gsl_vector* alloc_vector, size_t offset, size_t size) {
        _view = gsl_vector_subvector(alloc_vector, offset, size);
    }
    
    @property gsl_vector* vec() {
        return &_view.vector;
    }
}




void main() {
    auto vec = gsl_vector_alloc(1000);    

    auto vec1 = new SubVec(vec, 0, 10);
    auto vec2 = new SubVec(vec, 10, 10);

    gsl_vector_set_zero(vec1.vec);
    gsl_vector_set_zero(vec2.vec);
    gsl_vector_set(vec1.vec, 4, 152);
    gsl_vector_set(vec2.vec, 9, 144.222);

    foreach(i; 0 .. 10)
        writef("%s ", gsl_vector_get(vec1.vec, i));
    write("\n");

    foreach(i; 0 .. 10)
        writef("%s ", gsl_vector_get(vec2.vec, i));
    write("\n");
    writeln(vec1.vec().owner);
}
