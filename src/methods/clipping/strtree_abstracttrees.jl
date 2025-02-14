using SortTileRecursiveTree: STRtree, STRNode, STRLeafNode
@static if !Base.hasmethod(AbstractTrees.children, Tuple{STRtree})
    AbstractTrees.children(tree::STRtree) = AbstractTrees.children(tree.rootnode)
    AbstractTrees.parent(tree::STRtree) = nothing

    AbstractTrees.nodevalue(tree::STRtree) = AbstractTrees.nodevalue(tree.rootnode)

    # Implement the interface for general STRNodes

    AbstractTrees.children(node::STRNode) = node.children
    AbstractTrees.nodevalue(node::STRNode) = node.extent

    # Implement the interface for STRLeafNodes
    AbstractTrees.children(node::STRLeafNode) = STRLeafNode[]
    AbstractTrees.nodevalue(node::STRLeafNode) = reduce(Extents.union, node.extents)


    # Define the traits from AbstractTrees
    AbstractTrees.ParentLinks(::Type{<: STRtree}) = AbstractTrees.ImplicitParents()
    AbstractTrees.ParentLinks(::Type{<: STRNode}) = AbstractTrees.ImplicitParents()
    AbstractTrees.ParentLinks(::Type{<: STRLeafNode}) = AbstractTrees.ImplicitParents()

    AbstractTrees.SiblingLinks(::Type{<: STRtree}) = AbstractTrees.ImplicitSiblings()
    AbstractTrees.SiblingLinks(::Type{<: STRNode}) = AbstractTrees.ImplicitSiblings()
    AbstractTrees.SiblingLinks(::Type{<: STRLeafNode}) = AbstractTrees.ImplicitSiblings()

    AbstractTrees.ChildIndexing(::Type{<: STRtree}) = AbstractTrees.IndexedChildren()
    AbstractTrees.ChildIndexing(::Type{<: STRNode}) = AbstractTrees.IndexedChildren()
    # We don't define this trait for STRLeafNodes, because they have no children.

    # Type stability fixes

    AbstractTrees.NodeType(::Type{<:Union{STRNode, STRLeafNode, STRtree}}) = AbstractTrees.NodeTypeUnknown()
end