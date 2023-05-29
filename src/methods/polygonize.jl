
"""
    polygonize(A; minpoints=10)
    polygonize(xs, ys, A; minpoints=10)

Convert matrix `A` to polygons.

If `xs` and `ys` are passed in they are used as the pixel center points.

# Keywords
- `minpoints`: ignore polygons with less than `minpoints` points. 
"""
polygonize(A::AbstractMatrix; kw...) = polygonize(axes(A)..., A; kw...) 
function polygonize(xs, ys, A::AbstractMatrix; minpoints=10)
    contours = Iterators.map(contours(A)) do contour
            poly = map(contour) do xy
                x, y = Tuple(xy)
                Point2f(x + first(xs) - 1, y + first(ys) - 1)
            end
        end
    end 
    if minpoints > 1
        contours = Iterators.filter(contours) do contour
            length(contour) > minpoints
        end
       return map(Polygon, contours)
    else
        return map(Polygon, contours)
    end
end

# rotate direction clockwise
rot_clockwise(dir) = (dir) % 8 + 1
# rotate direction counterclockwise
rot_counterclockwise(dir) = (dir + 6) % 8 + 1

# move from current pixel to next in given direction
function move(pixel, image, dir, dir_delta)
    newp = pixel + dir_delta[dir]
    height, width = size(image)
    if (0 < newp[1] <= height) && (0 < newp[2] <= width)
        if image[newp] != 0
            return newp
        end
    end
    return CartesianIndex(0, 0)
end

# finds direction between two given pixels
function from_to(from, to, dir_delta)
    delta = to - from
    return findall(x -> x == delta, dir_delta)[1]
end

function detect_move(image, p0, p2, nbd, border, done, dir_delta)
    dir = from_to(p0, p2, dir_delta)
    moved = rot_clockwise(dir)
    p1 = CartesianIndex(0, 0)
    while moved != dir ## 3.1
        newp = move(p0, image, moved, dir_delta)
        if newp[1] != 0
            p1 = newp
            break
        end
        moved = rot_clockwise(moved)
    end

    if p1 == CartesianIndex(0, 0)
        return
    end

    p2 = p1 ## 3.2
    p3 = p0 ## 3.2
    done .= false
    while true
        dir = from_to(p3, p2, dir_delta)
        moved = rot_counterclockwise(dir)
        p4 = CartesianIndex(0, 0)
        done .= false
        while true ## 3.3
            p4 = move(p3, image, moved, dir_delta)
            if p4[1] != 0
                break
            end
            done[moved] = true
            moved = rot_counterclockwise(moved)
        end
        push!(border, p3) ## 3.4
        if p3[1] == size(image, 1) || done[3]
            image[p3] = -nbd
        elseif image[p3] == 1
            image[p3] = nbd
        end

        if (p4 == p0 && p3 == p1) ## 3.5
            break
        end
        p2 = p3
        p3 = p4
    end
end

"""
   contours(A::AbstractMatrix)

Returns contours as vectors of `CartesianIndex`.
"""
function contours(image::AbstractMatrix)
    nbd = 1
    lnbd = 1
    image = Float64.(image)
    contour_list = Vector{typeof(CartesianIndex[])}()
    done = [false, false, false, false, false, false, false, false]

    # Clockwise Moore neighborhood.
    dir_delta = (CartesianIndex(-1, 0), CartesianIndex(-1, 1), CartesianIndex(0, 1), CartesianIndex(1, 1), 
                 CartesianIndex(1, 0), CartesianIndex(1, -1), CartesianIndex(0, -1), CartesianIndex(-1, -1))

    height, width = size(image)

    for i = 1:height
        lnbd = 1
        for j = 1:width
            fji = image[i, j]
            is_outer = (image[i, j] == 1 && (j == 1 || image[i, j-1] == 0)) ## 1 (a)
            is_hole = (image[i, j] >= 1 && (j == width || image[i, j+1] == 0))

            if is_outer || is_hole
                # 2
                border = CartesianIndex[]
                from = CartesianIndex(i, j)

                if is_outer
                    nbd += 1
                    from -= CartesianIndex(0, 1)

                else
                    nbd += 1
                    if fji > 1
                        lnbd = fji
                    end
                    from += CartesianIndex(0, 1)
                end

                p0 = CartesianIndex(i, j)
                detect_move(image, p0, from, nbd, border, done, dir_delta) ## 3
                if isempty(border) ##TODO
                    push!(border, p0)
                    image[p0] = -nbd
                end
                push!(contour_list, border)
            end
            if fji != 0 && fji != 1
                lnbd = abs(fji)
            end

        end
    end

    return contour_list
end
