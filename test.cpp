class node {
public:
    int num;
    node* next;
    node(int val) {
        num = val;
        next = NULL;
    }
};

node* intersectionPresent(node* head1,node* head2) {
    node* d1 = head1;
    node* d2 = head2;

    while(d1 != d2) {
        d1 = d1 == NULL? head2:d1->next;
        d2 = d2 == NULL? head1:d2->next;
    }

    return d1;
}