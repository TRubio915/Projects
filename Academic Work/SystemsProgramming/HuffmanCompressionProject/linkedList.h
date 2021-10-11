typedef struct tnode
{
  int weight;
  char character;
  struct tnode* left;
  struct tnode* right;
  struct tnode* parent;
}huffTree;

typedef struct node
{
  huffTree* value;
  struct node* next;
} LinkedList;

LinkedList* llCreate();
int llIsEmpty(LinkedList* ll);
void llDisplay(LinkedList* ll);
void llAdd(LinkedList** ll, huffTree newValue);
void llFree(LinkedList* ll);
void listAddInOrder(LinkedList** ll, huffTree* newValue);
huffTree* removeFirst(LinkedList** ll);
