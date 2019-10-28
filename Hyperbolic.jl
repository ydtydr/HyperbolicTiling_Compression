#predefined functions
module Hyperbolic

    using StaticArrays
    using LinearAlgebra
    ##############################################################################
    # a point in S-dimensional hyperbolic space
    abstract type HyPoint{S} end;
    ##############################################################################
    # a point in the Lorentz model of hyperbolic space
    #
    #   S         the dimension of the hyperbolic space
    #   T         the underlying numeric type used for the model
    #   SS        S + 1
    struct HyPLorentz{S,T<:Number,SS} <: HyPoint{S}
        c :: SVector{SS,T};
        HyPLorentz(cc::SVector{SS,T}) where {SS,T<:Number} = new{SS-1,T,SS}(cc);
    end
    # distance between two points in hyperbolic space
    function dist(x::HyPLorentz{S,T,SS},y::HyPLorentz{S,T,SS}) where {S,T<:Number,SS}
        acc = x.c[1] * y.c[1];
        for i = 2:(S+1)
          acc -= x.c[i] * y.c[i];
        end
#         @assert(acc >= one(T));
        if acc < one(T)
            return zero(T)
        end
        return acosh(acc);
    end
    ##############################################################################
    # a point in the Poincare model of hyperbolic space
    #
    #   S         the dimension of the hyperbolic space
    #   T         the underlying numeric type used for the model
    struct HyPPoincare{S,T<:Number} <: HyPoint{S}
        c :: SVector{S,T};
        HyPPoincare(cc::SVector{S,T}) where {S,T<:Number} = new{S,T}(cc);
    end
    # distance between two points in hyperbolic space
    function dist(x::HyPPoincare{S,T}, y::HyPPoincare{S,T}) where {S,T<:Number}
        xnorm = 1 - dot(x.c,x.c);
        ynorm = 1 - dot(y.c,y.c);
#         @assert(xnorm > 0);
#         @assert(ynorm > 0);
        if xnorm<0 || ynorm<0
            return zero(T)
        end
        inn = 1 + 2 * dot(x.c-y.c,x.c-y.c) / (xnorm*ynorm);
        return acosh(inn);
    end
    ##############################################################################
    # a point in the L-tling model of hyperbolic space
    #
    #   F         the underlying floating point type used for the fundamental domain points
    #   I         the underlying big integer type used for the integer matrix
    #   c         the floating point in the fundamental domain
    #   g         the integer matrix
    struct HyPLTiling{F<:Number,I<:Integer} <: HyPoint{2}
        c :: SVector{3,F};
        g :: SMatrix{3,3,I};
        HyPLTiling(cc::SVector{3,F}, gg::SMatrix{3,3,I}) where {F<:Number,I<:Integer} = new{F,I}(cc,gg);
    end
    # distance between two points in hyperbolic space
    function dist(x::HyPLTiling{F,I},y::HyPLTiling{F,I}) where {F<:Number,I<:Integer}
        #Not implemented yet.
    end
    ##############################################################################
    # a point in the L-tling model of hyperbolic space with VBW encoding
    #
    #   F         the underlying floating point type used for the fundamental domain points
    #   c         the floating point in the fundamental domain
    #   vbw       the VBW encoding of the integer matrix
    struct HyPLTiling_VBW{F<:Number} <: HyPoint{2}
        c :: SVector{3,F};
        vbw :: BitArray{1}
        HyPLTiling_VBW(cc::SVector{3,F}, vbw_::BitArray{1}) where {F<:Number} = new{F}(cc, vbw_);
    end
    # distance between two points in hyperbolic space
    function dist(x::HyPLTiling_VBW{F},y::HyPLTiling_VBW{F}) where {F<:Number}
        #Not implemented yet.
    end
    ##############################################################################
    # convert from Lorentz model to Poincare model
    function HyPPoincare(x::HyPLorentz{S,T,SS}) where {S,T<:Number,SS}
        idxs = SVector{S}(2:SS);
        return HyPPoincare(x.c[idxs] / (1+x.c[1]));
    end
    # convert from Poincare model to Lorentz model
    function HyPLorentz(x::HyPPoincare{S,T}) where {S,T<:Number}
        normx2 = dot(x.c,x.c);
        @assert(normx2 < 1);
        return HyPLorentz(vcat(SVector((1 + normx2)/(1 - normx2)), x.c * (2 / (1 - normx2))));
    end
    # convert from L-Tiling model to Lorentz model
    function HyPLorentz(x::HyPLTiling{F,I}) where {F<:Number,I<:Integer}
        L = SMatrix{3,3,F}([sqrt(3*one(F)) 0 0;0 1 0;0 0 1]);
        LInv = SMatrix{3,3,F}([sqrt(3*one(F))/3 0 0;0 1 0;0 0 1]);
        return HyPLorentz(L*x.g*LInv*x.c);
    end
    # convert from L-Tiling_vbw model to Lorentz model
    function HyPLorentz(x::HyPLTiling_VBW{F}) where {F<:Number}
        L = SMatrix{3,3,F}([sqrt(3*one(F)) 0 0;0 1 0;0 0 1]);
        LInv = SMatrix{3,3,F}([sqrt(3*one(F))/3 0 0;0 1 0;0 0 1]);
        g = Re_Mat(x.vbw, BigInt)
        return HyPLorentz(L*g*LInv*x.c);
    end
    # convert from Lorentz model to L-Tiling model
    function HyPLTiling(x::HyPLorentz{2,T,3}, i::Type{I}) where {T<:Number, I<:Integer}
        L = SMatrix{3,3,T}([sqrt(3*one(T)) 0 0;0 1 0;0 0 1]);
        LInv = SMatrix{3,3,T}([sqrt(3*one(T))/3 0 0;0 1 0;0 0 1]);
        ga = SMatrix{3,3,I}([2 1 0;0 0 -1;3 2 0]);
        gb = SMatrix{3,3,I}([2 -1 0;0 0 -1;-3 2 0]);
        RVI = SMatrix{3,3,I}([1 0 0;0 1 0;0 0 1]);
        y = copy(x.c);
        y[1]=sqrt(1+y[2]^2+y[3]^2);
        numconut = 0;
        while 2*y[2]^2-y[3]^2-1>0 || 2*y[3]^2-y[2]^2-1>0
            numconut += 1;
            if y[2]<-abs(y[3])
                RVI = ga*RVI;
            elseif y[2]>abs(y[3])
                RVI = gb*RVI;
            elseif y[3]<-abs(y[2])
                RVI = gb^5*RVI;
            elseif y[3]>abs(y[2])
                RVI = ga^5*RVI;
            end
            y = L*RVI*R*x.c;
            y[1]=sqrt(1+y[2]^2+y[3]^2);
        end
        return HyPLTiling(y,RVI);
    end
    ##############################################################################
    #transform a string to bitarray
    function Str2BA(x::String)
        vbw = Int8[]
        for i=1:length(x)
            push!(vbw, parse(Int8,x[i]))
        end
        return BitArray(vbw)
    end
    #transform a bitarray to string
    function BA2Str(x::BitArray{1})
        str = ""
        for i=1:length(x)
            str *= string(Int(x[i]))
        end
        return str
    end
    ##############################################################################
    # test ga or gb type at the end of an orderstr
    #
    # orderstr, the VBW encoding string of an element in group
    function LastEle(orderstr::String)
        if parse(Int,orderstr[1])==0
            #test what's the generator order type at the left beginning, 0 for ga, 1 for gb
            abind = 0 #left beginning order is ga
        else
            abind = 1 #left beginning order is gb
        end
        nextind = 2 # current position in orderstr
        ord_len = 0 # the length of the generator order, according to its number
        while nextind<length(orderstr)
            if parse(Int,orderstr[nextind])==1
                #this order number is either 1 or 5, of length 2
                ord_len = 2
                nextind += ord_len
            else
                ord_len = 3#this order number is of length 3
                nextind += ord_len
            end
            abind = 1-abind # ga and gb appear alternatively
        end
        #output the last generator order at the end of orderstr, and length (to determine order number)
        return 1-abind, ord_len
    end
    ##############################################################################
    # Update the VBW encoding when an group element is multiplied with a new group element
    #
    # orderstr, the VBW encoding string of an element in group
    # whic, the matrix element to be multiplied to orderstr
    # 0 for ga^5, 1 for gb^5, 2 for gb, 3 for ga
    function Update_VBW(orderstr::String, whic::Int)
        @assert length(orderstr)>0
        abind, lastind = LastEle(orderstr)
        if lastind==2 && orderstr[end-1:end]=="10"
            #last element of orderstr is of order number 1
            if abind==0#0 for ga, last element of orderstr is ga
                if whic == 0#to be multiplied with ga^5
                    neworderstr = orderstr[1:end-2]
                elseif whic ==1#to be multiplied with gb^5
                    neworderstr = orderstr*"11"
                elseif whic ==2#to be multiplied with gb
                    neworderstr = orderstr*"10"
                elseif whic ==3#to be multiplied with ga
                    neworderstr = orderstr[1:end-2]*"001"
                end
            elseif abind==1#1 for gb
                if whic == 0#ga^5
                    neworderstr = orderstr*"11"
                elseif whic ==1#gb^5
                    neworderstr = orderstr[1:end-2]
                elseif whic ==2#gb
                    neworderstr = orderstr[1:end-2]*"001"
                elseif whic ==3#ga
                    neworderstr = orderstr*"10"
                end
            end
        elseif lastind==2 && orderstr[end-1:end]=="11"#order 5
            if abind==0#0 for ga
                if whic == 0#ga^5
                    neworderstr = orderstr[1:end-2]*"011"
                elseif whic ==1#gb^5
                    neworderstr = orderstr*"11"
                elseif whic ==2#gb
                    neworderstr = orderstr*"10"
                elseif whic ==3#ga
                    neworderstr = orderstr[1:end-2]
                end
            elseif abind==1#1 for gb
                if whic == 0#ga^5
                    neworderstr = orderstr * "11"
                elseif whic ==1#gb^5
                    neworderstr = orderstr[1:end-2]*"011"
                elseif whic ==2#gb
                    neworderstr = orderstr[1:end-2]
                elseif whic ==3#ga
                    neworderstr = orderstr* "10"
                end
            end
        elseif lastind==3 && orderstr[end-2:end]=="001"#order 2
            if abind==0#0 for ga
                if whic == 0#ga^5
                    neworderstr = orderstr[1:end-3]*"10"
                elseif whic ==1#gb^5
                    neworderstr = orderstr*"11"
                elseif whic ==2#gb
                    neworderstr = orderstr*"10"
                elseif whic ==3#ga
                    neworderstr = orderstr[1:end-3]*"010"
                end
            elseif abind==1#1 for gb
                if whic == 0#ga^5
                    neworderstr = orderstr * "11"
                elseif whic ==1#gb^5
                    neworderstr = orderstr[1:end-3]*"10"
                elseif whic ==2#gb
                    neworderstr = orderstr[1:end-3]*"010"
                elseif whic ==3#ga
                    neworderstr = orderstr* "10"
                end
            end
        elseif lastind==3 && orderstr[end-2:end]=="010"#order 3
            if abind==0#0 for ga
                if whic == 0#ga^5
                    neworderstr = orderstr[1:end-3]*"001"
                elseif whic ==1#gb^5
                    neworderstr = orderstr*"11"
                elseif whic ==2#gb
                    neworderstr = orderstr*"10"
                elseif whic ==3#ga
                    neworderstr = orderstr[1:end-3]*"011"
                end
            elseif abind==1#1 for gb
                if whic == 0#ga^5
                    neworderstr = orderstr * "11"
                elseif whic ==1#gb^5
                    neworderstr = orderstr[1:end-3]*"001"
                elseif whic ==2#gb
                    neworderstr = orderstr[1:end-3]*"011"
                elseif whic ==3#ga
                    neworderstr = orderstr* "10"
                end
            end
        elseif lastind==3 && orderstr[end-2:end]=="011"#order 4
            if abind==0#0 for ga
                if whic == 0#ga^5
                    neworderstr = orderstr[1:end-3]*"010"
                elseif whic ==1#gb^5
                    neworderstr = orderstr*"11"
                elseif whic ==2#gb
                    neworderstr = orderstr*"10"
                elseif whic ==3#ga
                    neworderstr = orderstr[1:end-3]*"11"
                end
            elseif abind==1#1 for gb
                if whic == 0#ga^5
                    neworderstr = orderstr * "11"
                elseif whic ==1#gb^5
                    neworderstr = orderstr[1:end-3]*"010"
                elseif whic ==2#gb
                    neworderstr = orderstr[1:end-3]*"11"
                elseif whic ==3#ga
                    neworderstr = orderstr* "10"
                end
            end
        end
        if length(neworderstr)==1
            return ""
        end
        return neworderstr
        #VBW encoding of orderstr multiplied with new matrix element
    end
    ##############################################################################
    # recover the matrix from the VBW encoding
    #
    # orderstr, the VBW encoding string of an element in group
    # i, to indicate the type of integer
    #
    function Re_Mat(orderstr::BitArray{1}, i::Type{I}) where {I<:Integer}
        orderstr = BA2Str(orderstr);
        RV= SMatrix{3,3,I}([1 0 0;0 1 0;0 0 1]);
        ga = SMatrix{3,3,I}([2 1 0;0 0 -1;3 2 0])
        gb = SMatrix{3,3,I}([2 -1 0;0 0 -1;-3 2 0])
        if length(orderstr) == 0
            return RV
        end
        if parse(Int,orderstr[1])==0
            abind = 0#0 for a, 1 for b
        else
            abind = 1#0 for a, 1 for b
        end
        nextind = 2
        while nextind<length(orderstr)
            if parse(Int,orderstr[nextind])==1
                if orderstr[nextind:nextind+1] == "10"
                    if abind == 0
                        RV = RV*ga
                    else
                        RV = RV*gb
                    end
                else
                    if abind == 0
                        RV = RV*ga^5
                    else
                        RV = RV*gb^5
                    end
                end
                nextind +=2
            else
                if orderstr[nextind:nextind+2] == "001"
                    if abind == 0
                        RV = RV*ga^2
                    else
                        RV = RV*gb^2
                    end
                elseif orderstr[nextind:nextind+2]=="010"
                    if abind == 0
                        RV = RV*ga^3
                    else
                        RV = RV*gb^3
                    end
                else
                    if abind == 0
                        RV = RV*ga^4
                    else
                        RV = RV*gb^4
                    end
                end
                nextind +=3
            end
            abind = 1-abind
        end
        return RV
    end
    ##############################################################################
    # compress a point in Lorentz model to L-tiling model, outputs the point in fundamental domain and the
    # matrix order VBW encoding, 1: "10",2: "001",3: "010",4: "011",5: "11"
    #
    #
    # x, a point in the Lorentz model
    # i, to indicate the type of integer
    # orderstr, the existent VBW encoding string of an element in group, can be empty or non-empty
    #
    #
    function Gen_VBW(x::HyPLorentz{2,T,3}, i::Type{I}, orderstr::String="") where {T<:Number, I<:Integer}
        L = SMatrix{3,3,T}([sqrt(3*one(T)) 0 0;0 1 0;0 0 1]);
        LInv = SMatrix{3,3,T}([sqrt(3*one(T))/3 0 0;0 1 0;0 0 1]);
        ga = SMatrix{3,3,I}([2 1 0;0 0 -1;3 2 0]);
        gb = SMatrix{3,3,I}([2 -1 0;0 0 -1;-3 2 0]);
        RVI = SMatrix{3,3,I}([1 0 0;0 1 0;0 0 1]);
        y = Array(copy(x.c));
        y[1]=sqrt(1+y[2]^2+y[3]^2);
        numconut = 0;
        while 2*y[2]^2-y[3]^2-1>0 || 2*y[3]^2-y[2]^2-1>0
            numconut += 1;
            if y[2]<-abs(y[3])
                if numconut==1
                    orderstr = orderstr*"011";
                else
                    orderstr = Update_VBW(orderstr,0);
                end
                RVI = ga*RVI;
            elseif y[2]>abs(y[3])
                if numconut==1
                    orderstr = orderstr*"111";
                else
                    orderstr = Update_VBW(orderstr,1);
                end
                RVI = gb*RVI;
            elseif y[3]<-abs(y[2])
                if numconut==1
                    orderstr = orderstr*"110";
                else
                    orderstr = Update_VBW(orderstr,2);
                end
                RVI = gb^5*RVI;
            elseif y[3]>abs(y[2])
                if numconut==1
                    orderstr = orderstr*"010";
                else
                    orderstr = Update_VBW(orderstr,3);
                end
                RVI = ga^5*RVI;
            end
            y = Array(L*RVI*LInv*x.c);
            y[1]=sqrt(1+y[2]^2+y[3]^2);
        end
        return HyPLTiling_VBW(SVector{3}(y), Str2BA(orderstr))#
    end
    ##########################
end
