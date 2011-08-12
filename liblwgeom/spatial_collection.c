#include "liblwgeom.h"
#include "spatial_collection.h"

SPATIAL_COLLECTION *
sc_create(COLLECTION_TYPE t,
		  PARAMETERS *params,
		  INCLUDES *inc,
		  EVALUATOR *eval)
{
	SPATIAL_COLLECTION *sc ;

	/* collection must have an includes */
	if (inc == NULL) return NULL ;

	sc = (SPATIAL_COLLECTION *)lwalloc(sizeof(SPATIAL_COLLECTION)) ;
	if ( sc != NULL ) {
		/* not a derived collection */
		sc->input1 = NULL ;
		sc->input2 = NULL ;

		/* copy provided pointers to the collection struct */
		sc->type = t ;
		sc->params = params ;
		sc->inclusion = inc ;
		sc->evaluator = eval ;

		/* set up reverse links */
		sc->inclusion->collection = sc ;
		if (sc->evaluator != NULL) {
			sc->evaluator->collection = sc ;
		}
	}

	return sc ;
}

SPATIAL_COLLECTION *
sc_twoinput_create(COLLECTION_TYPE t,
		           PARAMETERS *params,
		           INCLUDES   *inc,
		           EVALUATOR  *eval,
		           SPATIAL_COLLECTION *input1,
		           SPATIAL_COLLECTION *input2)
{
	SPATIAL_COLLECTION *sc ;

	sc = sc_create(t, params, inc, eval) ;
	if (sc != NULL) {
		sc->input1 = input1 ;
		sc->input2 = input2 ;
	}
	return sc ;
}

INCLUDES *
inc_create(PARAMETERS *params, INCLUDE_FN includes, INCLUDE_FN includesIndex)
{
	INCLUDES *inc ;

	/* includes must be specified */
	if (includes == NULL) return NULL ;
	inc = (INCLUDES *)lwalloc(sizeof(INCLUDES)) ;

	if (inc != NULL) {
		inc->params = params ;
		inc->includes = includes ;
		inc->includesIndex = includesIndex ;
		inc->collection = NULL ;
	}

	return inc ;
}

EVALUATOR *
eval_create(PARAMETERS *params,
		    EVALUATOR_FN evaluator,
		    EVALUATOR_FN evaluatorIndex,
		    int result_len)
{
	EVALUATOR *eval ;

	/* evaluator must be specified */
	if (evaluator == NULL) return NULL ;
	eval = (EVALUATOR *)lwalloc(sizeof(EVALUATOR)) ;

	if (eval != NULL) {
		eval->params = params ;
		eval->evaluator = evaluator ;
		eval->evaluatorIndex = evaluatorIndex ;
		eval->collection = NULL ;
		eval->result = val_create(result_len) ;
	}

	return eval ;
}

void
eval_destroy(EVALUATOR *eval)
{
	if (eval != NULL) {
		/* if there's an associated collection, un-associate it before freeing */
		if (eval->collection != NULL) {
			eval->collection->evaluator = NULL ;
		}
		lwfree(eval) ;
	}
}

void
inc_destroy(INCLUDES *inc)
{
	if (inc != NULL) {
		/* if there's an associated collection, un-associate it before freeing */
		if (inc->collection != NULL ) {
			inc->collection->inclusion = NULL ;
		}
		lwfree(inc) ;
	}
}

void
sc_destroy(SPATIAL_COLLECTION *sc)
{
	if (sc != NULL) {
		/* dis-associate the inclusion object */
		if (sc->inclusion != NULL) {
			sc->inclusion->collection = NULL ;
		}

		/* dis-associate the evaluation object */
		if (sc->evaluator != NULL) {
			sc->evaluator->collection = NULL ;
		}

		lwfree(sc) ;
	}
}

/* to compute the start of the data array */
#define PADDING (sizeof(double)-(sizeof(VALUE)%sizeof(double)))

VALUE *
val_create(int num_values)
{
	VALUE * val ;
	int i ;

	if (num_values <= 0) return NULL ;
	val = (VALUE *)lwalloc(sizeof(VALUE) + (sizeof(double)*num_values)+PADDING) ;

	/* initialize the data pointer to right after the struct */
	if (val != NULL) {
		val->data = (double *)(((char *)val)+sizeof(VALUE)+PADDING) ;
	}

	/* initialize data values */
	for (i=0; i<num_values; i++) {
		val->data[i] = 0.0 ;
	}

	val->length = num_values ;

	return val ;
}

void
val_destroy(VALUE *val)
{
	if (val != NULL) {
		lwfree(val) ;
	}
}

VALUE *
sc_evaluate(SPATIAL_COLLECTION *sc, LWPOINT *point)
{

	/* check for potential NULL values */
	if ((sc == NULL) || (point == NULL)) return NULL ;
	if (sc->evaluator == NULL) return NULL ;
	if (sc->evaluator->evaluate == NULL) return NULL ;

	return sc->evaluator->evaluate(sc->evaluator, point) ;
}

VALUE *
sc_evaluateIndex(SPATIAL_COLLECTION *sc, LWPOINT *point)
{

	/* check for potential NULL values */
	if ((sc == NULL) || (point == NULL)) return NULL ;
	if (sc->evaluator == NULL) return NULL ;
	if (sc->evaluator->evaluateIndex == NULL) return NULL ;

	return sc->evaluator->evaluateIndex(sc->evaluator, point) ;
}


int
sc_includes(SPATIAL_COLLECTION *sc, LWPOINT *point)
{
	/* check for potential NULL values */
	if ((sc == NULL) || (point == NULL)) return 0 ;
	if (sc->inclusion == NULL) return 0 ;
	if (sc->inclusion->includes == NULL) return 0 ;

	return sc->inclusion->includes(sc->inclusion, point) ;
}

int
sc_includesIndex(SPATIAL_COLLECTION *sc, LWPOINT *point)
{
	/* check for potential NULL values */
	if ((sc == NULL) || (point == NULL)) return 0 ;
	if (sc->inclusion == NULL) return 0 ;
	if (sc->inclusion->includesIndex == NULL) return 0 ;

	return sc->inclusion->includesIndex(sc->inclusion, point) ;
}

int
sc_hasValue(SPATIAL_COLLECTION *sc)
{
	if (sc == NULL) return 0 ;
	return (sc->type == SPATIAL_PLUS_VALUE) ;
}

int
sc_hasTwoInputs(SPATIAL_COLLECTION *sc)
{
	if (sc == NULL) return 0 ;
	return (sc->input1 != NULL) && (sc->input2 != NULL) ;
}
