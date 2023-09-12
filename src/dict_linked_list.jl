using DataStructures

"""
Data structure containing a sorted linked list and a dictionary pointing to the
nodes in the list. The types K and V are the key types for the dictionary and
the value types for the linked list respectively. Assume the type V is totally
ordered.
"""
mutable struct DictLinkedList{K, V}
    # Sorted linked list
    list :: MutableLinkedList{V}
    # Dictionary with with values pointing nodes in the linked list
    dict :: Dict{K, DataStructures.ListNode{V}}
    # total order on V
    comp :: Function

    DictLinkedList{K, V}(comp) where {K, V} = new(
        MutableLinkedList{V}(),
        Dict{K, DataStructures.ListNode{V}}(),
        comp)
end


"""
Make a DictLinkedList iterable by iterating over the underlying linked list
"""
function Base.iterate(dll :: DictLinkedList{K, V}, state = dll.list.node) where {K, V}
    root_node = dll.list.node
    state.next === root_node ? nothing : (state.next.data, state.next)
end


"""
Insert a value into a DictLinkedList, preserving the order of the underlying
list and adding a key to the underlying dictionary.
"""
function insert!(dll :: DictLinkedList{K, V}, key :: K, value :: V) :: Nothing where {K, V}
    # check the key does not already exist in the dictionary
    haskey(dll.dict, key) && error("Duplicate key")

    # Pointer to the element immediately before where value should be inserted
    elem_before = dll.list.node
    if dll.comp(elem_before.prev.data, value)
        # can insert at the end of the list
        elem_before = elem_before.prev
    else
        next = iterate(dll)
        while next !== nothing
            (x, current) = next
            if dll.comp(x, value)
                elem_before = current
                next = iterate(dll, current)
            else
                break
            end
        end
    end
    # insert the value here
    new_node = DataStructures.ListNode{V}(value)
    new_node.prev = elem_before
    new_node.next = elem_before.next
    elem_before.next.prev = new_node
    elem_before.next = new_node
    dll.list.len += 1

    # add an entry to the dictionary
    dll.dict[key] = new_node
    return nothing
end


"""
Remove a key from the underlying dictionary of a DictLinkedList and remove the
corresponding value from the underlying linked list
"""
function remove!(dll :: DictLinkedList{K, V}, key :: K) :: Nothing where {K, V}
    !haskey(dll.dict, key) && error("Key does not exist")
    node = dll.dict[key]
    delete!(dll.dict, key)
    node.prev.next = node.next
    node.next.prev = node.prev
    dll.list.len -= 1
    return nothing
end

"""
Implementation of Base.empty function for DictLinkedList. We will simply check
that the underlying linked list is empty
"""
Base.empty(dll :: DictLinkedList{K, V}) where {K, V} = iterate(dll) === nothing


Base.length(dll :: DictLinkedList{K, V}) where {K, V} = dll.list.len


==(dll1 :: DictLinkedList{K, V}, dll2 :: DictLinkedList{K, V}) where {K, V} = (
    dll1.list == dll2.list)


function filter!(dll :: DictLinkedList{K, V}, f :: Function) :: Nothing where {K, V}
    remove_stack = Stack{K}(dll.list.len)
    for (k, v) in dll.dict
        f(v.data) ? continue : push!(remove_stack, k)
    end
    for k in remove_stack
        remove!(dll, k)
    end
end


"""
Check if a key is in the underlying dictionary of a DictLinkedList
"""
Base.contains(dll :: DictLinkedList{K, V}, key :: K) where {K, V} = haskey(dll.dict, key)


"""
Find the value corresponding to the given key
"""
lookup(dll :: DictLinkedList{K, V}, key :: K) where {K, V} = dll.dict[key].data


"""
Replace the value associated with a given key. Assume that the order of the list
is preserved.
"""
function replace!(dll :: DictLinkedList{K, V}, key :: K, value :: V) where {K, V}
    # @assert dll.comp(value, lookup(dll, key)) && dll.comp(lookup(dll, key), value)
    dll.dict[key].data = value
end


"""
Get the last element of the underlying list
"""
last(dll :: DictLinkedList{K, V}) where {K, V} = dll.list.node.prev.data


"""
Implementation of the Base.map function for DictLinedList
"""
function Base.map(f :: Function, dll :: DictLinkedList{K, V}) where {K, V}
    dll_ = deepcopy(dll)
    next = iterate(dll_)
    while next !== nothing
        (x, current) = next
        current.data = f(x)
        next = iterate(dll_, current)
    end
    return dll_
end


"""
Show a DictLinkedList
"""
function show(io :: IO, dll :: DictLinkedList{K, V}) where {K, V}
    println(io, "DictLinkedList{$K, $V}")
    println(io, '\t', dll.dict)
    println(io, '\t', dll.list)
end
