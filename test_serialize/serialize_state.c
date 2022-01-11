//-----------------------------------------
// Notes: fbuffer.h is not included in:
//        /usr/local/include/msgpack.h
//-----------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>     // for strchr()
#include <sys/stat.h>   // to get filesize
#include <msgpack.h>
#include "msgpack/fbuffer.h"        // See note above
#include "../include/bmi.h"
#include "../include/serialize_state.h"

//------------------------------------------------
// Now computing these from filesize of ser_file
//------------------------------------------------
// static const int BUFFER_SIZE = 65536;
// static const int UNPACKED_BUFFER_SIZE = 131072;
// static const int BUFFER_SIZE = 1024;
// static const int UNPACKED_BUFFER_SIZE = 2048;

/*
These functions are meant to provide a "framework utility"
that can serialize, save to file, and then deserialize all
of the state variables of a model written in C.

The function "serialize()" uses new BMI functions to retrieve
the model's state variables called:  get_state_var_ptrs(),
get_state_var_names(), get_state_var_types(), and
get_state_var_lengths().

The function "deserialize_to_state()" deserializes the saved model
(from a file) and then uses the new BMI function set_state_var()
to set all state variable values into a new instance of the model.
This new instance could be running on another node.

The function "compare_states()" compares all of the state
variables of two instances of the given model, to confirm that
they are identical after calling "deserialize_to_state()".

These new functions could be used for the problem of
load balancing or for recovery after hardware failures.

See this documentation for msgpack:
https://github.com/msgpack/msgpack-c/wiki/v2_0_c_overview

See this msgpack example:
https://blog.gypsyengineer.com/en/security/msgpack-fuzzing.html
*/
//-----------------------------------------------------------------------
// Functions in this file:
//     get_file_size()
//     get_state_var_count()
//     get_state_var_names()
//     get_state_var_types()
//     get_state_var_lengths()
//     get_state_var_ptrs()
//    
//     serialize()
//     deserialize_to_state()
//     compare_states()
//-----------------------------------------------------------------------
int get_file_size(const char *ser_file, unsigned long int *file_size){

    //----------------------------------------
    // Get the file size, to set BUFFER_SIZE
    //----------------------------------------
    struct stat st;
    stat(ser_file, &st);
    *file_size = st.st_size;
    return 0;
}

//-----------------------------------------------------------------------
int get_state_var_count(Bmi *model1, int *count){

    if (!model1){
        return BMI_FAILURE;   
    }
    
    char role[] = "all";
    model1->get_model_var_count( model1, role, count);
    
    //puts("In get_state_var_count()...");    
    //printf("   count = %d\n", *count);
    //puts("");
    return BMI_SUCCESS;    
}
//-----------------------------------------------------------------------
int get_state_var_names(Bmi *model1, char **names){

    // Should we add "struct' before "Bmi" above?
    if (!model1){
        return BMI_FAILURE;   
    }

    //------------------------
    // Do this in the caller
    //------------------------
    //char **names = NULL;
    //names = (char**) malloc (sizeof(char *) * n_state_vars);
    //for (int i=0; i<n_state_vars; i++){
    //    names[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);

    
    char role[] = "all";
    model1->get_model_var_names( model1, role, names);

    //--------------
    // For testing
    //--------------
    //size_t n_state_vars = sizeof(names) / sizeof(names[0]); 
    ////int n_state_vars = sizeof(names) / sizeof(names[0]);
           
    //puts("In get_state_var_names()...");
    //printf("   role     = %s\n", role);
    //printf("   names[0] = %s\n", names[0]);
    //printf("   names[1] = %s\n", names[1]);
    //printf("   names[2] = %s\n", names[2]);
    
    //printf("   sizeof(names)    = %zud\n", sizeof(*names) );
    //printf("   sizeof(names[0]) = %zud\n", sizeof(*names[0]) );
    //printf("   n_state_vars = %zud\n", n_state_vars);
    // printf("   n_state_vars = %d\n", n_state_vars);
    //puts("");
    
    return BMI_SUCCESS;    
}

//-----------------------------------------------------------------------
int get_state_var_types(Bmi *model1, char **names, char **types){

    if (!model1){
        return BMI_FAILURE;
    }

    // char *type;  //#### Doesn't work now.
    char type[BMI_MAX_TYPE_NAME];   // This works.
    int n_state_vars;
    get_state_var_count(model1, &n_state_vars);
    
    //puts("In get_state_var_types()...");
    //printf("n_state_vars = %d\n", n_state_vars);
    
    //int n_state_vars = sizeof(names) / sizeof(names[0]);
        
    // Added "names" as argument.  #################
    //char role[] = "all";
    //model1->get_model_var_names( model1, role, names, &n_state_vars);
    
    for (int i = 0; i < n_state_vars; i++) {
        model1->get_var_type( model1, names[i], type );
        strncpy( types[i], type, BMI_MAX_TYPE_NAME);
        // types[i] = type;  //#### Doesn't work now.
        //printf("types[i] = %s\n", types[i]);
    }
  
    return BMI_SUCCESS;
}

//-----------------------------------------------------------------------
int get_state_var_lengths(Bmi *model1, char **names, unsigned int lengths[]){

    if (!model1){
        return BMI_FAILURE;   
    }

    int length, n_state_vars, itemsize, nbytes;
    get_state_var_count(model1, &n_state_vars);
    
    for (int i = 0; i < n_state_vars; i++) {
        //----------------------------------------    
        // Method that uses new get_var_length()
        //----------------------------------------
        model1->get_var_length( model1, names[i], &length );
        lengths[i] = length;
        
        //----------------------------------
        // Method without get_var_length()
        //----------------------------------
        //model1->get_var_itemsize( model1, names[i], &itemsize);
        //model1->get_var_nbytes(   model1, names[i], &nbytes);
        //lengths[i] = (nbytes / itemsize);

        //printf("itemsize = %d\n", itemsize);
        //printf("nbytes   = %d\n", nbytes);
        //printf("length   = %d\n", lengths[i]);
        //printf("\n");
    }

    //puts("In get_state_var_lengths()...");
    //printf("   names[0]   = %s\n", names[0]);
    //printf("   lengths[0] = %d\n", lengths[0]);
    //puts("");   
    return BMI_SUCCESS;
}

//-----------------------------------------------------------------------
int get_state_var_ptrs(Bmi *model1, char **names, void *ptrs[]){

    if (!model1){
        return BMI_FAILURE;   
    }

    void *ptr;
    int n_state_vars;
    get_state_var_count(model1, &n_state_vars);
    
    // int n_state_vars = sizeof(names) / sizeof(names[0]);

    // Added "names" as argument. #################
    //char role[] = "all";
    //model1->get_model_var_names( model1, role, names, &n_state_vars);
 
    //puts("In get_state_var_ptrs()...");   
    //printf("   n_state_vars = %d\n", n_state_vars);
    //printf("   names[0] = %s\n", names[0]);
    
    for (int i = 0; i < n_state_vars; i++) {
        //printf("   name = %s\n", names[i]);
        model1->get_value_ptr( model1, names[i], &ptr );
        ptrs[i] = ptr;
    }
  
    return BMI_SUCCESS;
}

//---------------------------------------------------------
// We don't need this now.  Can just use get_file_size().
//---------------------------------------------------------
// int get_buffer_size(Bmi* model1, unsigned int buffer_size){
//

//     int n_state_vars;
//     model1->get_state_var_count(model1, &n_state_vars);
// 
//     unsigned int length, lengths[ n_state_vars ], type_size;
//     unsigned int unpacked_buffer_size = 0;
//     int verbose = 1;
//     int i;
//
//     char type[BMI_MAX_TYPE_NAME];
//     char **types = NULL;
//     types = (char**) malloc (sizeof(char *) * n_state_vars);  
//     for (i=0; i<n_state_vars; i++){
//         types[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
//     }
// 
//     if (verbose){ puts("Calling BMI.get_state_var_types()..."); }
//     model1->get_state_var_types(model1, types);
// 
//     if (verbose){ puts("Calling BMI.get_state_var_lengths()..."); }
//     model1->get_state_var_lengths(model1, lengths); 
// 
//     for (i=0; i<n_state_vars; i++){
//         type   = types[i];
//         length = lengths[i];
//  
// 		if (strcmp(type, "int") == 0){
// 			type_size = sizeof( int );
// 		} else if (strcmp(type, "float") == 0){
// 			type_size = sizeof( float );
// 		} else if (strcmp(type, "double") == 0){
// 			type_size = sizeof( double );
// 		} else if (strcmp(type, "string") == 0){
// 			type_size = sizeof( char ); 
// 		} else if (strcmp(type, "FILE") == 0){
// 			type_size = sizeof( FILE );
// 		} else{
// 			printf("  WARNING: Unknown type = %s\n", types[i] );
// 			type_size = 4;
// 		}
//                    
//         unpacked_buffer_size += (size * type_size);
//     }
//    
//     buffer_size = unpacked_buffer_size;
//     
//     //--------------------------------   
//     // Free memory for string arrays
//     //--------------------------------
//     if (verbose){ puts("Freeing memory..."); }
//     for (i=0; i<n_state_vars; i++){
//         free (types[i]);
//     }
//     free (types);
// 
//     if (verbose){ puts("Finished computing buffer size."); puts(""); }
//     return 0;
//                  
// } 

//-----------------------------------------------------------------------
int serialize(Bmi* model1, const char *ser_file) {

    int n_state_vars, result;
    result = get_state_var_count( model1, &n_state_vars);
    printf("In serialize()...\n");
    printf("   n_state_vars = %d.\n", n_state_vars);
    
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_count(); returning.");
        return BMI_FAILURE;
    }
    // model1->get_state_var_count(model1, &n_state_vars);  //##############
    
    FILE *fp = fopen(ser_file, "w+");
    int i, j, n_bytes;
    void *ptr_list[ n_state_vars ];
    unsigned int length, lengths[ n_state_vars ];
    int verbose = 1;
    int    i_val;
    long   li_val;
    float  f_val;
    double d_val;
  
    char **names = NULL;
    names = (char**) malloc (sizeof(char *) * n_state_vars);
    for (i=0; i<n_state_vars; i++){
        names[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }
    
    // char *type;  // Doesn't work now.  #######
    char type[BMI_MAX_TYPE_NAME];   // This works.
    char **types = NULL;
    types = (char**) malloc (sizeof(char *) * n_state_vars);  
    for (i=0; i<n_state_vars; i++){
        types[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }

    //------------------------------------------------------
    // Get required info about the model's state variables
    //------------------------------------------------------ 
    if (verbose){ puts("Calling get_state_var_names() in serialize()..."); }
    result = get_state_var_names(model1, names);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_names(); returning.");
        return BMI_FAILURE;
    }
    //if (verbose){ puts("Calling BMI.get_state_var_names()..."); }
    //model1->get_state_var_names(model1, names);
    
    if (verbose){ puts("Calling get_state_var_types() in serialize()..."); }
    result = get_state_var_types(model1, names, types);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_types(); returning.");
        return BMI_FAILURE;
    }
    //if (verbose){ puts("Calling BMI.get_state_var_types()..."); }
    //model1->get_state_var_types(model1, types);

    if (verbose){ puts("Calling get_state_var_lengths() in serialize()..."); }
    result = get_state_var_lengths(model1, names, lengths);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_lengths(); returning.");
        return BMI_FAILURE;
    }
    //if (verbose){ puts("Calling BMI.get_state_var_lengths()..."); }
    //model1->get_state_var_lengths(model1, lengths);       

    if (verbose){ puts("Calling get_state_var_ptrs() in serialize()..."); }
    result = get_state_var_ptrs(model1, names, ptr_list);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_ptrs(); returning.");
        return BMI_FAILURE;
    }
    //if (verbose){ puts("Calling BMI.get_state_var_ptrs()..."); }
    //model1->get_state_var_ptrs(model1, ptr_list);

    //--------------------------------------------
    // Prepare to write serialized state to file
    // msgpack_fbuffer_write needs fbuffer.h
    //--------------------------------------------
    msgpack_packer pk;
    msgpack_packer_init(&pk, fp, msgpack_fbuffer_write);

    //----------------------------------------------
    // Prepare to write serialized state to buffer
    //----------------------------------------------
    // msgpack_sbuffer* buffer = msgpack_sbuffer_new();
    // msgpack_packer* pk = msgpack_packer_new(buffer, msgpack_sbuffer_write);
  
    //------------------------------------------------------
    // Note: strings are serialized with 2 commands, like:
    //------------------------------------------------------
    // msgpack_pack_str(&pk, 7);
    // msgpack_pack_str_body(&pk, "example", 7);

    //--------------------------------------
    // Note: booleans are serialized like:
    //--------------------------------------
    // msgpack_pack_true(&pk);
    // msgpack_pack_false(&pk);

    //----------------------------------------------
    // Dereference the pointers & serialize values
    // Note that we aren't using the names here.
    //---------------------------------------------- 
    if (verbose){
        puts("Serializing the model state variables...");
        printf("  Number of state vars = %d\n", n_state_vars);
    } 
    for (i=0; i<n_state_vars; i++){
        length = lengths[i];
        strncpy( type, types[i], BMI_MAX_TYPE_NAME);
        // type = types[i];  //#### Doesn't work now.
        // printf("  types[i] = %s\n", types[i]);
        printf("  i = %d, type = %s, length = %u\n", i, type, length );

        if (ptr_list[i] == NULL){
            if (verbose){
                printf("  NULL pointer encountered for i = %d.\n", i); }
        } 

        //-------------------------------
        // Does type name contain "*" ?
        //-------------------------------
        // printf("BEFORE type = %s\n", type);
        if (strchr(type, '*') != NULL){
            type[strlen(type)-1] = '\0';  }   // remove last char
        // printf("AFTER type  = %s\n", type);
        // printf("\n");

        //---------------------------------------------
        // Is this state var a single value or array?
        //---------------------------------------------
        if (length == 1){
            if (strcmp(type, "int") == 0){
                msgpack_pack_int(&pk, *(int *)ptr_list[i]);
            } else if (strcmp(type, "long") == 0){
                msgpack_pack_long(&pk, *(long *)ptr_list[i]);
            } else if (strcmp(type, "float") == 0){
                msgpack_pack_float(&pk, *(float *)ptr_list[i]);  
            } else if (strcmp(type, "double") == 0){
                msgpack_pack_double(&pk, *(double *)ptr_list[i]);
            } else if (strcmp(type, "string") == 0){
                // Note:  Need ptr_list[i] without * here.
                msgpack_pack_str_body(&pk, ptr_list[i], length);
            } else if (strcmp(type, "char") == 0){
                msgpack_pack_str_body(&pk, ptr_list[i], length);      
            } else if (strcmp(type, "FILE") == 0){
                // nil is an object pointer to nothing
                msgpack_pack_nil(&pk);  // Need something; will this work?
            } else{
                printf("  WARNING: Unknown type = %s\n", types[i] );
                msgpack_pack_nil(&pk);  // Need something; will this work?
            }
        } else{
            //---------------------------------------
            // This state var is an array or string
            //---------------------------------------
            if (strcmp(type, "string") == 0){
                msgpack_pack_str(&pk, length);
            } else if (strcmp(type, "char") == 0){
                msgpack_pack_str(&pk, length);
            } else if (strcmp(type, "FILE") != 0){
                msgpack_pack_array(&pk, length);
                // n_bytes = length * sizeof( double );
                // msgpack_pack_bin(&pk, n_bytes);
            }

            //-----------------------------------------------------
            // NOTE! Typecast ptr first, then add offset,
            //       into array, then dereference the ptr.
            //       CORRECT:    *( ((double *)ptr) + j)
            //       INCORRECT:  *( (double *)(ptr + j))
            //       INCORRECT:  *( (double *)ptr + j)  ??
            //       INCORRECT:  *(double *)ptr + j 
            //---------------------------------------------
            if (strcmp(type, "int") == 0){
                for (j=0; j<length; j++){
                    i_val = *( ((int *)ptr_list[i]) + j);
                        msgpack_pack_int(&pk, i_val); }
            } else if (strcmp(type, "long") == 0){
                for (j=0; j<length; j++){
                    li_val = *( ((long *)ptr_list[i]) + j);
                    msgpack_pack_long(&pk, li_val); }    //#####################
            } else if (strcmp(type, "float") == 0){
                for (j=0; j<length; j++){
                    f_val = *( ((float *)ptr_list[i]) + j);
                    msgpack_pack_float(&pk, f_val); }
            } else if (strcmp(type, "double") == 0){
                for (j=0; j<length; j++){
                    d_val = *( ((double *)ptr_list[i]) + j);
                    msgpack_pack_double(&pk, d_val); }
            } else if (strcmp(type, "string") == 0){
                // Note:  Need ptr_list[i] without * here.
                msgpack_pack_str_body(&pk, ptr_list[i], length);
            } else if (strcmp(type, "char") == 0){
                // Note:  Need ptr_list[i] without * here.
                msgpack_pack_str_body(&pk, ptr_list[i], length);     
            } else if (strcmp(type, "FILE") == 0){
                for (j=0; j<length; j++){
                    // nil is an object pointer to nothing
                     msgpack_pack_nil(&pk); } // Need something; will this work?
            } else{
                printf("  WARNING: Unknown type = %s\n", types[i] );
                msgpack_pack_nil(&pk);  // Need something; will this work?
            }
        }
    }
        
    //---------------------------------
    // If writing to file, close file
    //---------------------------------
    if (verbose){ puts("Closing output file..."); }
    fclose(fp);

    //------------------------------------  
    // If writing to buffer, free memory
    //------------------------------------
    // msgpack_sbuffer_free(buffer);
    // msgpack_packer_free(pk);
  
    //--------------------------------   
    // Free memory for string arrays
    //--------------------------------
    if (verbose){ puts("Freeing memory..."); }
    for (i=0; i<n_state_vars; i++){
        free (names[i]);
        free (types[i]);
    }
    free (names);
    free (types);

    if (verbose){ puts("Finished serializing."); puts(""); }
    return 0;
}
//--------------------------------------------------------------------
int deserialize_to_state(const char *ser_file, Bmi* model2, int print_obj) {

    //-----------------------------------------------------------
    // This online reference shows how to unpack the serialized
    // buffer one object at a time, as done here.
    //-----------------------------------------------------------
    // https://github.com/msgpack/msgpack-c/blob/c_master/
    //         example/user_buffer_unpack.c

    //-----------------------------------------------------------
    // This online reference shows how to save the serialized
    // buffer to a file and read it back in, also done here.
    //-----------------------------------------------------------
    // https://blog.gypsyengineer.com/en/security/msgpack-fuzzing.html
    
    //------------------------
    // Additional references
    //--------------------------------------------------------------------       
    // https://stackoverflow.com/questions/15393838/
    //         how-do-i-unpack-and-extract-data-properly-using-msgpack-c
    // https://stackoverflow.com/questions/12431441/messagepack-c-api/
    //         12581029#12581029
    //--------------------------------------------------------------------
    int verbose   = 1;
    // int print_obj = 1;
    int     i_val,  *i_arr, j, result;
    short   si_val, *si_arr;
    long    li_val, *li_arr;
    float   f_val,  *f_arr;
    double  d_val,  *d_arr;
    void    *ptr;

    //------------------------------------     
    // Get all the total state var count
    //------------------------------------
    int n_state_vars;
    result = get_state_var_count(model2, &n_state_vars);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_count(); returning.");
        return BMI_FAILURE;
    }
      
    //------------------------------   
    // Get all the state var names
    //------------------------------
    // char *name;
    char name[BMI_MAX_VAR_NAME];
    char **names = NULL;
    names = (char**) malloc (sizeof(char *) * n_state_vars);
    for (int j=0; j<n_state_vars; j++){
        names[j] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }
    result = get_state_var_names(model2, names);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_names(); returning.");
        return BMI_FAILURE;
    }
    
    //------------------------------     
    // Get all the state var types
    //------------------------------
    // char *type;
    char type[BMI_MAX_TYPE_NAME];
          
    char *sval;
    char **types = NULL;
    types = (char**) malloc (sizeof(char *) * n_state_vars);  
    for (int j=0; j<n_state_vars; j++){
        types[j] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }
    result = get_state_var_types(model2, names, types);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_types(); returning.");
        return BMI_FAILURE;
    }

    //--------------------------------     
    // Get all the state var lengths
    //--------------------------------
    unsigned int length, lengths[ n_state_vars ];
    result = get_state_var_lengths(model2, names, lengths);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_lengths(); returning.");
        return BMI_FAILURE;
    }
    
    //-------------------------------------
    // Get the file size; set buffer_size
    //-------------------------------------
    unsigned long int file_size, buffer_size, unpacked_buffer_size;
    get_file_size( ser_file, &file_size );
    buffer_size = file_size;
    unpacked_buffer_size = 2 * buffer_size;
    char inbuffer[buffer_size];
    char unpacked_buffer[unpacked_buffer_size]; 
//     if (verbose){
//         printf("Buffer_size = %lu\n", buffer_size);
//         printf("Unpacked buffer_size = %lu\n", unpacked_buffer_size);
//     }  
    // char inbuffer[BUFFER_SIZE];
    // char unpacked_buffer[UNPACKED_BUFFER_SIZE];
        
    FILE *fp = fopen(ser_file, "rb");
    int i = 0;
    size_t off = 0;
    size_t len = 0;

    msgpack_unpacked unpacked;
    msgpack_unpack_return ret;
    msgpack_unpacked_init(&unpacked);
    if (verbose){ puts("Deserializing source model state vars..."); } 
    
    //------------------------------------------------------------
    // In online ref 1, buffer is passed as an argument (buf).
    // In online ref 2, buffer is read from ser_file (inbuffer).
    //------------------------------------------------------------
    // In online ref 1, len is passed as an argument.
    // In online ref 2, len=read is return value from fread.
    //------------------------------------------------------------
    // Online ref 1: "result" = Online ref 2: "unpacked".
    // Online ref 1: "buf" = Online ref 2: "inbuffer"  ????
    // Online ref 1: "len" = Online ref 2: "read"
    //------------------------------------------------------------
    // len = fread(inbuffer, sizeof(char), BUFFER_SIZE, fp); 
    len = fread(inbuffer, sizeof(char), buffer_size, fp); 
    ret = msgpack_unpack_next(&unpacked, inbuffer, len, &off);

    //-------------------------------------------
    // Deserialize and set state vars in model2
    //-------------------------------------------
    while (ret == MSGPACK_UNPACK_SUCCESS) {
        msgpack_object obj = unpacked.data;

        //#### type = types[i];
        //#### name = names[i];
        strncpy( type, types[i], BMI_MAX_TYPE_NAME);
        strncpy( name, names[i], BMI_MAX_VAR_NAME);
        length = lengths[i];

        // printf("type = %s\n", type);
        if (strchr(type, '*') != NULL){
            type[strlen(type)-1] = '\0';  // remove last char
        }

        if (length == 1){
            //------------------------------------------------
            // Note:  These commented lines don't work, even
            //        though we typecast in set_state_var().
            //------------------------------------------------
            // ptr = &obj;
            // model2->set_state_var(model2, ptr, i );
            //--------------------------------------------                  
            if (strcmp(type, "int") == 0){
                i_val = (int)obj.via.i64;
                model2->set_value(model2, name, &i_val );
            } else if (strcmp(type, "short") == 0){        
                si_val = (short)obj.via.i64;
                model2->set_value(model2, name, &si_val );
            } else if (strcmp(type, "long") == 0){        
                li_val = (long)obj.via.i64;
                model2->set_value(model2, name, &li_val );
            } else if (strcmp(type, "float") == 0){        
                f_val = (float)obj.via.f64;
                model2->set_value(model2, name, &f_val );
            } else if (strcmp(type, "double") == 0){   
                d_val = (double)obj.via.f64;
                model2->set_value(model2, name, &d_val );
            } else if (strcmp(type, "char") == 0){
                model2->set_value(model2, name, &obj );
            } else if (strcmp(type, "string") == 0){
                model2->set_value(model2, name, &obj );
            } else{
                 model2->set_value(model2, name, &obj );
            }
        } else{
            //--------------------------------------------------
            // Copy values from msgpack array object into an
            // array of the correct type, then pass pointer to
            // that array to set_state_var(). (This works.)
            //--------------------------------------------------
            // Couldn't figure out a clean way to access the
            // msgpack array elements with pointer math. This
            // seems to be because it allows mixed types.
            // Options:  via.u64, via.i64, via.f64, via.str.
            //--------------------------------------------------
            if (strcmp(type, "string") == 0){
                // Works for title & subcat strings in TOPMODEL.
                ptr = obj.via.array.ptr;
                model2->set_value(model2, name, ptr );
            //--------------------------------------------------
            } else if (strcmp(type, "char") == 0){
                ptr = obj.via.array.ptr;
                model2->set_value(model2, name, ptr );
            //--------------------------------------------------            
            } else if (strcmp(type, "int") == 0){
                i_arr = (int*) malloc(length * sizeof( int ));
                for (int j=0; j<length; j++){
                    i_arr[j] = obj.via.array.ptr[j].via.i64; }
                model2->set_value(model2, name, i_arr );
                free(i_arr);
            //--------------------------------------------------
            } else if (strcmp(type, "long") == 0){
                li_arr = (long*) malloc(length * sizeof( long ));   //##########
                for (int j=0; j<length; j++){
                    li_arr[j] = obj.via.array.ptr[j].via.i64; }
                model2->set_value(model2, name, li_arr );
                free(li_arr);
            //--------------------------------------------------
            } else if (strcmp(type, "short") == 0){
                si_arr = (short*) malloc(length * sizeof( short ));
                for (int j=0; j<length; j++){
                    si_arr[j] = obj.via.array.ptr[j].via.i64; }
                model2->set_value(model2, name, si_arr );
                free(si_arr);
            //--------------------------------------------------
            } else if (strcmp(type, "float") == 0){
                f_arr = (float*) malloc(length * sizeof( float ));
                for (int j=0; j<length; j++){
                    f_arr[j] = obj.via.array.ptr[j].via.f64; }
                model2->set_value(model2, name, f_arr );
                free(f_arr);
            //--------------------------------------------------
            } else if (strcmp(type, "double") == 0){
                d_arr = (double*) malloc(length * sizeof( double ));
                for (int j=0; j<length; j++){
                    d_arr[j] = obj.via.array.ptr[j].via.f64; }
                model2->set_value(model2, name, d_arr );
                free(d_arr);
            }
            
            // For testing
            //printf("### obj.via.array.ptr[0].via.f64 = %f\n", obj.via.array.ptr[0].via.f64);
            //printf("### obj.via.array.ptr[1].via.f64 = %f\n", obj.via.array.ptr[1].via.f64);
            //printf("### obj.via.array.ptr[2].via.f64 = %f\n", obj.via.array.ptr[2].via.f64);
            //printf("\n");
        }
              
        if (print_obj){
            printf("Object no %d:\n", i);
            msgpack_object_print(stdout, obj);
            printf("\n");
              
            //------------------------------
            // This just prints obj again
            //------------------------------
            // msgpack_object_print_buffer(unpacked_buffer, UNPACKED_BUFFER_SIZE, obj);
            // printf("%s\n", unpacked_buffer);        
        }       
        //-----------------------------------------------------
        i++;
        ret = msgpack_unpack_next(&unpacked, inbuffer, len, &off);
    }   // end of while loop

    //--------------------------------------------------------
    // Did we unpack the expected number of state variables?
    //--------------------------------------------------------
    if (i < n_state_vars){
        printf("WARNING: Expected %d state variables \n", n_state_vars);
        printf("         But unpacked only %d vars. \n", i);
        printf("         BUFFER_SIZE may be too small.\n");
        printf("");
    } else{
        printf("Unpacked %d state variables.\n", i);
        printf("");     
    }

    //if (ret == MSGPACK_UNPACK_CONTINUE) {
    //    printf("Every msgpack_object in the buffer was processed.\n");
    //}
    //else if (ret == MSGPACK_UNPACK_PARSE_ERROR) {
    //    printf("The data in the buffer has an invalid format.\n");
    //}

    msgpack_unpacked_destroy(&unpacked);
    fclose(fp);  // Try moving this up to just after fread.

    //--------------------------------   
    // Free memory for string arrays
    //--------------------------------
    if (verbose){ puts("Freeing memory...");}
    for (i=0; i<n_state_vars; i++){
        free (names[i]);
        free (types[i]);
    }
    free (names);
    free (types);
    
    if (verbose){ puts("Finished deserializing."); puts("");}
    return 0;
}

//------------------------------------------------------------------------
int compare_states(Bmi* model1, Bmi* model2){

    int    verbose = 1;
    int    i, j, match, n_state_vars, result;
    int    err_count = 0;
    int    i_val1, i_val2;
    short  si_val1, si_val2;
    long   li_val1, li_val2;
    float  f_val1, f_val2;
    double d_val1, d_val2;
    char   c_val1, c_val2;

    //--------------------------------------
    // Get total number of state variables
    //--------------------------------------
    result = get_state_var_count(model1, &n_state_vars);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_count(); returning.");
        return BMI_FAILURE;
    }

    //----------------------------------------
    // Get the state variable internal names
    //----------------------------------------
    // char *name;
    char name[BMI_MAX_VAR_NAME];
    char **names = NULL;
    names = (char**) malloc (sizeof(char *) * n_state_vars);
    for (i=0; i<n_state_vars; i++){
        names[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }
    result = get_state_var_names(model1, names);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_names(); returning.");
        return BMI_FAILURE;
    }
      
    //------------------------------------  
    // Get the state variable data types
    //------------------------------------
    // char *type;
    char type[BMI_MAX_TYPE_NAME];
    char **types = NULL;
    types = (char**) malloc (sizeof(char *) * n_state_vars);  
    for (i=0; i<n_state_vars; i++){
        types[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
    }
    result = get_state_var_types(model1, names, types);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_types(); returning.");
        return BMI_FAILURE;
    }
    
    //---------------------------------  
    // Get the state variable lengths
    //---------------------------------
    unsigned int length, lengths[ n_state_vars ];
    result = get_state_var_lengths(model1, names, lengths);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_lengths(); returning.");
        return BMI_FAILURE;
    }

    //-------------------------------------        
    // Get pointers to Model 1 state vars
    //-------------------------------------
    void *ptr_list1[ n_state_vars ];
    result = get_state_var_ptrs(model1, names, ptr_list1);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_ptrs(); returning.");
        return BMI_FAILURE;
    }

    //-------------------------------------
    // Get pointers to Model 2 state vars
    //-------------------------------------  
    void *ptr_list2[ n_state_vars ];
    result = get_state_var_ptrs(model2, names, ptr_list2);
    if (result == BMI_FAILURE){
        puts("ERROR in get_state_var_ptrs(); returning.");
        return BMI_FAILURE;
    }

    //--------------------------------
    // Loop over all state variables
    //--------------------------------
    for (int i=0; i<n_state_vars; i++){
        match = 1;
        // type = types[i];
        // name  = names[i];  //########## This works, but "type" doesn't?
        strncpy( type, types[i], BMI_MAX_TYPE_NAME);
        strncpy( name, names[i], BMI_MAX_VAR_NAME);
        length = lengths[i];
        // printf("i = %d\n", i);  //##############
      
        //-------------------------------
        // Does type name contain "*" ?
        //-------------------------------
        if (strchr(type, '*') != NULL){
            type[strlen(type)-1] = '\0';  }   // remove last char
          
        //----------------------------------------------
        // NOTE! Typecast ptr first, then add offset,
        //       into array, then dereference the ptr.
        //       CORRECT:    *( ((double *)ptr) + j)
        //       INCORRECT:  *( (double *)(ptr + j))
        //       INCORRECT:  *( (double *)ptr + j)  ??
        //       INCORRECT:  *(double *)ptr + j 
        //---------------------------------------------
        if (strcmp(type, "int") == 0){
            for (j=0; j<length; j++){
                i_val1 = *( ((int *)ptr_list1[i]) + j);
                i_val2 = *( ((int *)ptr_list2[i]) + j);
                if (i_val1 != i_val2){
                    printf("Mismatch: i = %d, j = %d\n", i, j);
                    printf("i_val1 = %d, i_val2 = %d\n", i_val1, i_val2);
                    printf("var_name = %s\n\n", name);
                    match = 0;}
            }
        } else if (strcmp(type, "short") == 0){      //################
            for (j=0; j<length; j++){
                si_val1 = *( ((short *)ptr_list1[i]) + j);
                si_val2 = *( ((short *)ptr_list2[i]) + j);
                if (si_val1 != si_val2){
                    printf("Mismatch: i = %d, j = %d\n", i, j);
                    printf("si_val1 = %hd, si_val2 = %hd\n", si_val1, si_val2);
                    printf("var_name = %s\n\n", name);
                    match = 0; }
            }
        } else if (strcmp(type, "long") == 0){      //################
            for (j=0; j<length; j++){
                li_val1 = *( ((long *)ptr_list1[i]) + j);
                li_val2 = *( ((long *)ptr_list2[i]) + j);
                if (li_val1 != li_val2){
                    printf("Mismatch: i = %d, j = %d\n", i, j);
                    printf("li_val1 = %ld, li_val2 = %ld\n", li_val1, li_val2);
                    printf("var_name = %s\n\n", name);
                    match = 0; }
            }
        } else if (strcmp(type, "float") == 0){
            for (j=0; j<length; j++){
                f_val1 = *( ((float *)ptr_list1[i]) + j);
                f_val2 = *( ((float *)ptr_list2[i]) + j);
                if (f_val1 != f_val2){
                    printf("Mismatch: i = %d, j = %d\n", i, j);
                    printf("f_val1 = %f, f_val2 = %f\n", f_val1, f_val2);
                    printf("var_name = %s\n\n", name);
                    match = 0; }
            }
        } else if (strcmp(type, "double") == 0){
            for (j=0; j<length; j++){
                d_val1 = *( ((double *)ptr_list1[i]) + j);
                d_val2 = *( ((double *)ptr_list2[i]) + j);
                if (d_val1 != d_val2){
                    printf("Mismatch: i = %d, j = %d\n", i, j);
                    printf("d_val1 = %f, d_val2 = %f\n", d_val1, d_val2);
                    printf("var_name = %s\n\n", name);
                    match = 0; }
            }
        } else if (strcmp(type, "string") == 0){
              if (strcmp(ptr_list1[i], ptr_list2[i]) != 0){
                    printf("Mismatch: i = %d\n", i);
                    printf("str1 = %s\n", ptr_list1[i]);
                    printf("str2 = %s\n", ptr_list2[i]);                          
                    printf("var_name = %s\n\n", name);
                    match = 0;
              }
        } else if (strcmp(type, "char") == 0){
              if (strcmp(ptr_list1[i], ptr_list2[i]) != 0){
                    printf("Mismatch: i = %d\n", i);
                    printf("str1 = %s\n", ptr_list1[i]);
                    printf("str2 = %s\n", ptr_list2[i]);                          
                    printf("var_name = %s\n\n", name);
                    match = 0;
              }
        } else if (strcmp(type, "FILE") == 0){
              j = 0;  // do nothing 
        } else{
              printf("  Unknown type = %s\n", type );
        }

        //---------------------------------      
        // Did the state variables match?
        //---------------------------------
        if (match == 0){
            err_count++;
        } else{
            printf("Match: i = %d, name = %s\n", i, name) ;        
        }
    }

    //--------------------------------   
    // Free memory for string arrays
    //--------------------------------
    if (verbose){ puts("\nFreeing memory...");}
    for (i=0; i<n_state_vars; i++){
        free (names[i]);
        free (types[i]);
    }
    free (names);
    free (types);
  
    printf("Number of state variable differences between\n");
    printf("   the two model states = %d\n", err_count);
    printf("Finished comparing model states.\n\n");

    return 0;
}

//------------------------------------------------------------------------


